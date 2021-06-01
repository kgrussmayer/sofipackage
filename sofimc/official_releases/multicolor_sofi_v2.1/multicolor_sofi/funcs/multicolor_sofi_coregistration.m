function [imstack1crop, imstack2crop, calibration, rect, settings] = ...
    multicolor_sofi_coregistration(imstack1, imstack2, settings)
% Coregister two stacks of images. imstack2 is transformed on imstack1. 
% Boths image stacks are cropped to cover the same field of view.
%
% Inputs:
% imstack1      stack of images from first camera (rows, columns, frames)
% imstack2      stack of images from second camera (rows, columsn, frames)
% settings      [struct] settings for processing steps
% 
% Outputs:
% imstack1crop  coregistered and cropped image stack 1 
%               (rows, columns, frames)
% imstack2crop  coregistered and cropped image stack 2
%               (rows, columns, frames)
% calibration   [struct] tranformation used for coregistration
% rect          corregistration mask
% settings      [struct] settings for processings steps

% Copyright © 2018-2020 Tomas Lukes, Adrien Descloux
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

disp('Compute stacks STD')
st1 = std(single(imstack1),[],3);
st2 = std(single(imstack2),[],3);

%% Register channels based on calibration file with beads
if settings.cal.beads == 1  
    %%% Load the calibration files
    imc1 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc1, ...
        settings.io.fileExtension]);
    imc2 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc2, ...
        settings.io.fileExtension,]);
    if strcmp(settings.io.dataType, 'exper')
        imc2 = flip(imc2,2);
    end

    %%% Run calibration (register images from 2 color channels)
    disp('Compute calibration')
    settings.cal.roix = []; settings.cal.roiy = [];
    calibration = multicolor_sofi_calibration(...
        imc1,imc2,settings.sys,settings.cal,settings.io); 
end
  
%% Cross-correlations based channel registration
if settings.cal.beads == 0
    BORDER_SIZE = 50;
    disp('Compute correct coregistration')
    t1 = st1(BORDER_SIZE:end-BORDER_SIZE,BORDER_SIZE:end-BORDER_SIZE); 
    t2 = st2(BORDER_SIZE:end-BORDER_SIZE,BORDER_SIZE:end-BORDER_SIZE);
    % suppress salt&pepper noise
    t1 = medfilt2(t1,[2 2]); t2 = medfilt2(t2,[2 2]); 
    % suppress background
    t1 = t1-imgaussfilt(single(t1),5); t2 = t2-imgaussfilt(single(t2),5);

    temp = xcorr2(t1,t2);
    temp((size(temp,1)+1)/2,:) = 0;
    temp(:,(size(temp,2)+1)/2) = 0;
    [~,ind] = max(temp(:));
    [y,x] = ind2sub(size(temp),ind);
    shiftx = -(x-(size(temp,2)+1)/2);
    shifty = -(y-(size(temp,1)+1)/2);
    calibration = [];
end

%% Transform stack of images with the chosen calibration method
disp('Transform of imstack2')
Nframes = size(imstack1,3);
imstack2t = zeros(size(imstack1),'single');
% apply transform to data
if settings.cal.beads == 0 % using cross-correlation of raw data
    for ii = 1:Nframes
        im_mov = imstack2(:,:,ii);    
        im1t = imtranslate(im_mov,[-shiftx -shifty]);
        imstack2t(:,:,ii) = single(im1t);
    end
    settings.cal.shiftx = -shiftx;
    settings.cal.shifty = -shifty; % 
else % using beads calibration
    for ii = 1:Nframes
        im_mov = imstack2(:,:,ii);   
        im1t = imwarp(im_mov,calibration.tf{1},'OutputView',imref2d(size(im_mov)));
        imstack2t(:,:,ii) = single(im1t);
    end
    settings.cal.shiftx = -calibration.tf{1}.A(1);
    settings.cal.shifty = -calibration.tf{1}.B(1);
end

% compute coregistration mask and proper crop
mask = not(im1t(:,:,1) == 0);
[r,c] = find(mask>0);
rect = [min(c)+1 min(r)+1 max(c)-min(c)-2 max(r)-min(r)-2];
% get even dimensions for further SOFI processing 
if mod(rect(3),2) == 0; rect(3) = rect(3)-1; end
if mod(rect(4),2) == 0; rect(4) = rect(4)-1; end

settings.cal.roix = [rect(1) rect(3)];
settings.cal.roiy = [rect(2) rect(4)];

% crop image stack 
imstack1crop = imstack1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);
% crop image stack 
imstack2crop = imstack2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);

