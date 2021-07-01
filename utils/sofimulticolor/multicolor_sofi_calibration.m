function cal = multicolor_sofi_calibration(img_fix, img_mov, sys, cal, io)
% Run calibration. Calculate transformation matrices in between img_fix and
% img_mov using calibration stack of beads.
%
% Inputs:
% img_fix       stack of images from first camera (rows, columns, frames)
% img_mov       stack of images from second camera (rows, columsn, frames)
% sys           [struct] system settings
% cal           [struct] calibration settings
% io            [struct] input/output settings
%
% Outputs:
% calibration   [struct] tranformation in between img_fix and img_mov

% Copyright © 2018-2019 Tomas Lukes
% École Polytechnique Fédérale de Lausanne,
% Laboratory of Nanoscale Biology, http://lben.epfl.ch/
% tomas.lukes@epfl.ch

% This file is part of multicolorSOFI.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if cal.figs ==1
    % output path
    outputPath = [io.pnc,filesep,'testAlignment',filesep,io.fnc1];
    
    if (~exist(outputPath,'dir'))
        mkdir(outputPath);
    end
end

% load and prepare calibration data
% use only part of the images as specified in cal.roi
if ~isempty(cal.roiy)
    img_fix = img_fix(cal.roiy,:,:);
    img_mov = img_mov(cal.roiy,:,:);
end

if ~isempty(cal.roix)
    img_fix = img_fix(:,cal.roix,:);
    img_mov = img_mov(:,cal.roix,:);
end

% estimate transformation, test regitration

[~,indcom]=max(squeeze(mean(max(img_fix))).*squeeze(mean(max(img_mov))));
% median filter in a 3x3 neighbourhood
im_fix=medfilt2(img_fix(:,:,indcom),[3 3]);
im_mov=medfilt2(img_mov(:,:,indcom),[3 3]);

im_fix = max(im_fix,[],3);
im_mov = max(im_mov,[],3);

% extracting the center of gravities of the beads
out=struct;
out.ru=5;
sys.bg=1;
out.nh=0;

out=hriSegmentation(double(im_fix),cal.bgth,cal.logsize,out);
out=hriFilterSegments(double(im_fix),cal.aupl,cal.alol,sys,out);
cog_fix=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];

out=hriSegmentation(double(im_mov),cal.bgth,cal.logsize,out);
out=hriFilterSegments(double(im_mov),cal.aupl,cal.alol,sys,out);
cog_mov=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];

% identify corresponding beads

% spatial cross-correlation algorithm to determine shift of coordinates
[my, mx] = ccrShiftEstimation(im_fix,im_mov,10); % ! switched mx and my

cog_common_fix=zeros(1,2);
cog_common_mov=zeros(1,2);

%identify corresponding center of gravities

for k=1:size(cog_fix,1)
    d=(cog_mov-repmat(cog_fix(k,:)-[mx my],[size(cog_mov,1) 1]));
    d=sqrt(d(:,1).^2+d(:,2).^2);
    [p,ind]=min(d);
    if(p<cal.px_tol)
        cog_common_fix(end+1,:)=cog_fix(k,:);
        cog_common_mov(end+1,:)=cog_mov(ind,:);
    end
end

cog_common_fix(1,:)=[];  %first element is irrelevant (0,0)
cog_common_mov(1,:)=[];

if cal.figs ==1
    figure(),
    imshowpair(log(double(im_fix)+100), log(double(im_mov)+100),'Scaling','joint'); %green: first image, magenta: second image
    title('detected beads, log scale')
    hold on;
    scatter(cog_fix(:,2),cog_fix(:,1),'x','MarkerEdgeColor',[0 1 0]);
    scatter(cog_mov(:,2),cog_mov(:,1),'x','MarkerEdgeColor',[1 0 1]);
    scatter(cog_common_fix(:,2),cog_common_fix(:,1),'o','MarkerEdgeColor',[0 1 0]);
    scatter(cog_common_mov(:,2),cog_common_mov(:,1),'o','MarkerEdgeColor',[1 0 1]);
    legend('fix cog','mov cog','fix cog pair','mov cog pair');
    % save the current figure in tif file
    fh = gcf;saveFigure(fh,outputPath,'detectedBeads','tif');
end

% switch x,y coordinates
ccm=fliplr(cog_common_mov);
ccf=fliplr(cog_common_fix);

tform = fitgeotrans(ccf,ccm,'polynomial', 2);
[x,y] = transformPointsInverse(tform, cog_common_mov(:,2),cog_common_mov(:,1));
d=([y x]-cog_common_fix);
disp(['avg reg. error using polynomial geom. transform: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);

for ii = 1:length(sys.orders)
    tf{ii} = fitgeotrans(ii*ccm,ii*ccf,'polynomial', 2);
end

if cal.figs ==1
    moving_reg = imwarp(im_mov, tf{1},'OutputView',imref2d(size(im_fix)));
    figure(),
    imshowpair(log(double(im_fix)+100), log(double(moving_reg)+100),'Scaling','joint'); %green: first image, magenta: second image
    title('registered planes, log scale')
    % save the current figure in tif file
    fh = gcf;saveFigure(fh,outputPath,'registeredPlanes','tif');
    save([outputPath,'calib.mat'],'tf');
end

cal.tf = tf;
% eof