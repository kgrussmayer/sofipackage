function [sofi]=sofiLinearize3Dpad(sofi,dec)
% Deconvolve flat cumulants and linearize the brightness
% by taking the order-th root of the deconvolved cumulants. The result
% is reconvolved with a realistic point-spread function.
% Tomas Lukes, tomas.lukes@epfl.ch

%Inputs:
% sofi   flattened cumulants 
% dec    struct
%       .fwhm      PSF fwhm 
%       .noise     lative noise level                         
%       .orders    cumulant orders                    
%       .iter      maximum number of iterations for deconvoluion   
%
% Output:
% sofi  3D cumulant image

dims=numel(dec.fwhm);
fwhm=dec.fwhm/sqrt(8*log(2));

for order=dec.orders
    sigma=fwhm*sqrt(order);
    PSF=gaussian_psf(sigma(1));
    PSF=PSF(:)*gaussian_psf(sigma(2));
    images=abs(sofi{order});

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

   if strcmp(dec.psfmodel,'airy')
       if isfield(dec,'resx') == 1 && dec.useEstimate == 1
            PSF = airyPSF3D([dec.resx(order) dec.resx(order) ...
                fwhm(3)*sqrt(8*log(2))*(order)]).^order;
       else
            PSF = airyPSF3D(fwhm*sqrt(8*log(2))*(order)).^order;
       end
   end

    for k=1:k
        tic
        image=single(images(:,:,:,k));
        scale=max(image(:));
        if dec.axialMirror == 1
            temp = image;
            image = zeros([size(image,1),size(image,2),size(image,3)+6*order]);
            image(:,:,(1+3*order):end-3*order) = temp;
            image(:,:,1:3*order) = temp(:,:,3*order:-1:1);
            image(:,:,end-3*order+1:end) = temp(:,:,end:-1:end-3*order+1);
        end
        try
            % try to run the deconvolution on gpu 
            %(requires modification of deconvlucy.m function)
            image=(scale*gather(deconvlucy(gpuArray(image/scale),gpuArray(single(PSF)),dec.iter))).^(1/order);
        catch
            % if calling deconvlucy with gpuArray arguments failed, run on cpu
            image=scale*gather(deconvlucy(image/scale,PSF,dec.iter)).^(1/order);   
        end
        gpuDevice(1); % reset gpu, might avoid later cuda error 
        
        if dec.reconvolve == 1
            images(:,:,:,k)=convn(image,PSF,'same');
        else
             if dec.axialMirror == 1
                 images(:,:,:,k)=image(:,:,(1+3*order):end-3*order);
             else
                images(:,:,:,k)=image;
             end
        end
        t = toc;
        disp(['Order: ',num2str(order),', stack: ',num2str(k),', processing time: ',num2str(t),'(sec)'])
    end
    
    sofi{order}=reshape(images,sizes);
end
% eof