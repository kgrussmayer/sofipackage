function [umx_c2, im_dec, im_lin, mean1, mean2, mc2] = multicolor_sofi_interface(...
    imstack1, imstack2, settings)
% Run multicolor SOFI processing steps: registration, cumulant calculation,
% spectral unmixing, flattening, deconvolution, linearisation.
%
% Inputs:
% imstack1      stack of images from first camera (rows, columns, frames)
% imstack2      stack of images from second camera (rows, columsn, frames)
% settings      [struct] settings for processing steps
%
% Outputs:
% umx_c2        unmixed second order-cross cumulant images - sofi2 images
%               (rows, cols, frames, color channels)
% im_dec        deconvolved sofi images (rows, cols, frames, color channels)
% im_lin        linearized sofi images (rows, cols, frames, color channels)
% mean1         temporal average over frames of the corregistered stack1
% mean2         temporal average over frames of the corregistered stack2
% mc2           second order-cross cumulant images before unmixing - sofi2 images
%               (rows, cols, frames, color channels)

% Copyright © 2018-2020 Tomas Lukes
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

% Coregister image stacks
[imstack1crop, imstack2crop, calibration, rect, settings] = ...
    multicolor_sofi_coregistration(imstack1, imstack2, settings);
Nframes = size(imstack1,3);

% calculate the mean transformed image
mean1 = mean(double(imstack1crop), 3);
mean2 = mean(double(imstack2crop), 3);

% Calculate MC SOFI cumulants
% cropped raw data (1st stack) and transformed and cropped data (2nd stack)
disp('Compute SOFI cumulants')
orders = settings.sys.orders;
subseqlength = settings.sys.subseqlength;
if isempty(settings.sys.start), settings.sys.start = 1; end
start = settings.sys.start;

Nss=floor((Nframes-start+1)/subseqlength); % number of subsequences

for ns=1:Nss
    disp(['processing Nss:',num2str(ns)]);
    fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
    sofi=mcSofiCumulants(permute(cat(4,imstack1crop(:,:,fr),...
        imstack2crop(:,:,fr)),[1 2 4 3]),[],[],[],orders);
    if ns==1
        c=sofi;
        c=cellfun(@(x)repmat(x,[1 1 1 Nss]),c,'UniformOutput',0);
    else
        for io=orders
            c{io}(:,:,:,ns)=sofi{io};
        end
    end
end

% Calculate 2D cross-cumulants of non-transformed raw data (2nd stack)
disp('Compute SOFI cumulant of imstack2')
for ns=1:Nss
    disp(['processing Nss:',num2str(ns)]);
    fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
    sofi=mcSofiCumulantsSingle(imstack2(:,:,fr),[],[],[],orders);
    if ns==1
        c2=sofi;
        c2=cellfun(@(x)repmat(x,[1 1 Nss]),c2,'UniformOutput',0);
    else
        for io=orders
            c2{io}(:,:,ns)=sofi{io};
        end
    end
end

% Flattening before unmixing if chosen
if strcmp(settings.sys.flattening, 'before')
    disp('Cumulant flattening')
    c = sofiAllFlatten3D(c, orders);
    
    % flattening of 2D cross-cumulants of non-transformed raw data
    c2 = sofiAllFlatten3D(c2,orders);
end

% Transform second channel cross-cumulants and combine cumulants
disp('Transform and crop second SOFI channel')
c2t=c2;
off = [1 4 6];
for io=orders
    if settings.cal.beads == 0
        % transform and crop
        c2t{io} = imtranslate(c2{io},io.*[settings.cal.shiftx settings.cal.shifty]); % -cal.tf{1}.B(1)
        c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-off(io):io*rect(2) + ...
            io*rect(4)-ceil(3*io/2)-off(io),...
            io*rect(1)+floor(3*io/2)-off(io):io*rect(1)+io*rect(3)-ceil(3*io/2)-off(io),:);
    else% using beads calibration
        % transform and crop
        c2t{io} = imwarp(c2{io},calibration.tf{io},'OutputView',imref2d(size(c2{io})));
        c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-off(io):io*rect(2) + ...
            io*rect(4)-ceil(3*io/2)-off(io),...
            io*rect(1)+floor(3*io/2)-off(io):io*rect(1)+io*rect(3)-ceil(3*io/2)-off(io),:);
    end
end

disp('Combine cumulants')
tmp=cellfun(@(x)reshape(x(:,:,:),size(x,1),size(x,2),1,[]),c2t,'UniformOutput',0);
for io=orders
    c{io}(:,:,end,:)=tmp{io};
end

% Spectral unmixing - using 2nd order cumulant - 3 colors
disp('Spectral unmixing')
% load the transmission coefficient T (proportion in transmission channel)
T1 = settings.mc.T1;
T2 = settings.mc.T2;
T3 = settings.mc.T3;

R = [T1^2 T2^2 T3^2;
    T1*(1-T1) T2*(1-T2) T3*(1-T3);
    (1-T1)^2 (1-T2)^2 (1-T3)^2];
Rinv=R^(-1);

mc2=c{2};
umx1c2=squeeze(Rinv(1,1)*mc2(:,:,3,:)+Rinv(1,2)*mc2(:,:,2,:)+Rinv(1,3)*mc2(:,:,1,:));
umx2c2=squeeze(Rinv(2,1)*mc2(:,:,3,:)+Rinv(2,2)*mc2(:,:,2,:)+Rinv(2,3)*mc2(:,:,1,:));
umx3c2=squeeze(Rinv(3,1)*mc2(:,:,3,:)+Rinv(3,2)*mc2(:,:,2,:)+Rinv(3,3)*mc2(:,:,1,:));

% Flattening after unmixing if chosen
if strcmp(settings.sys.flattening, 'after')
    % flattening after unmixing
    umx_c2 = {umx1c2, umx2c2, umx3c2};
    umx_c2 = sofiAllFlattenMC(umx_c2,2);
    umx_c2 = cat(4, umx_c2{1}, umx_c2{2}, umx_c2{3});
end

% Deconvolution and Linearization of second order unmixed and flattened SOFI images
disp('Deconvolution')
order=2;
im_lin = [];
% psfs with zero padding for different deconvolutions
h = fspecial('Gaussian',29,settings.dec.fwhm(1)*sqrt(order)/sqrt(8*log(2)));

for n=1:size(umx1c2,3)
    for index_colors = 1:size(umx_c2,4)
        tmp = umx_c2(:,:,n,index_colors);
        tmp(tmp<0) = 0;
        scale=max(tmp(:));
        tmp = tmp./scale;
        if strcmp(settings.deconv_method,'lucy')
            im_dec = deconvlucy(tmp,h,settings.dec.iter);
            
        elseif strcmp(settings.deconv_method,'bsofi')
            % use always 50 iterations as in the bSOFI paper
            % (over-deconvolve and reconvole in the next step)
            im_dec = deconvlucy(tmp,h,50);
            
        else
            disp('Deconvolution method not implemented')
        end
        % optional linearization
        if settings.dec.lin == 1
            im_lin(:,:,n,index_colors) = scale*(im_dec.^settings.sys.lincoeff(2));
        else
            im_lin(:,:,n,index_colors) = scale*(im_dec);
        end
        
        if strcmp(settings.deconv_method,'bsofi')
            % reconvolve the image
            psf=gaussian_psf(settings.dec.fwhm/sqrt(8*log(2))); %
            psf=psf(:)*gaussian_psf(settings.dec.fwhm/sqrt(8*log(2)));
            im_lin(:,:,n,index_colors)=convn(im_lin(:,:,n,index_colors),psf,'same');
        end
    end
end
% eof