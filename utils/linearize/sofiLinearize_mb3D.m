function [sofi]=sofiLinearize_mb3D(sofi,dec)
% Deconvolve and denoise flat cumulants using Bregman iteration based 
% SOFI adapted deconvolution. Linearize the brightness of the deconvolved 
% cumulants.
% Tomas Lukes, tomas.lukes@epfl.ch
%
%Inputs:
% sofi      Flat cumulants (see sofiFlatten)
% dec   structure
% ..fwhm      PSF diameter (see sofiFlatten)
% noise     Noise level                         {1%}
% .orders    Cumulant orders                     {all}
% .iter      Maximum deconvolution iterations    {10}
%Output:
% sofi      Deconvolved linearized cumulant images


dims=numel(dec.fwhm);
fwhm=dec.fwhm/sqrt(8*log(2));

for order=dec.orders
    sigma=fwhm;
    PSF=gaussian_psf(sigma(1));
    PSF=PSF(:)*gaussian_psf(sigma(2));
    images=sofi{order};
    
    sizes=size(images);
    if dims == 3
        [x,y]=size(PSF);
        PSF=PSF(:)*gaussian_psf(sigma(3));
        PSF=reshape(PSF,x,y,size(PSF,2));
        k=prod(sizes(4:end));
    else
        k=prod(sizes(3:end));
        images=reshape(images,[sizes(1:2) 1 k]);
    end
    for k=1:k
        tic
        disp(k)
        img=images(:,:,:,k);
        scale=max(img(:));
        
        img_mproj = max(img,[],3);
        lambda = estimate_noise(img_mproj./max(img_mproj(:)));
        if lambda < 0.001; lambda = 0.001; end
        if lambda > 0.1; lambda = 0.1; end
        
        [sy, sx, sz] = size(img);
        
        dec.coef_A = 1;
        dec.lambda = 0.02;
        dec.apodize = 1;
        Aaray = [dec.iter,dec.apodize,dec.lambda,dec.fcpx,dec.coef_A,1];
        
        masksign = ones(sy,sx,sz);
        masksign (img<0)=-1;
        img = abs(img);
        img = (img./scale).^dec.lincoeffs(order);
        img = masksign.*img; % preserve sign
        
        for x=1:order
            for y=1:order
                temp=img(y:order:sizes(1),x:order:sizes(2),:);
                temp = SOFI_X7(temp,PSF,Aaray);
                img(y:order:sizes(1),x:order:sizes(2),:) = temp;
            end
        end
        img = abs(img);
        img = scale*img;
        
        images(:,:,:,k)=img;
        toc
    end
    sofi{order}=reshape(images,sizes);
end
% eof