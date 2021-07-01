function [sofi]=sofiLinearize(sofi,fwhm,orders,iter,settings,lincoeffs)
% Deconvolve flat cumulants and linearize the brightness
% by taking the order-th root of the deconvolved cumulants. The result
% is reconvolved with a realistic point-spread function.
% Tomas Lukes, tomas.lukes@epfl.ch
%
%Inputs:
% sofi          Flat cumulants (see sofiFlatten)
% fwhm          PSF diameter (see sofiFlatten)
% orders        Cumulant orders                     {all}
% iter          Maximum deconvolution iterations    {10}
% lincoeffs     Linearization coefficients
%
%Output:
% sofi          Linarized cumulant images

if nargin < 3 || isempty(orders)
   orders=find(~cellfun(@isempty,sofi));
end

if isempty(orders)
   return;
end

if nargin < 6
   lincoeffs = 1./orders;
end

if nargin < 4 || isempty(iter)
   iter=10;
end

if strcmp(settings.dec.psfmodel,'airy')
    psf = airyPSF(fwhm(1));
else % gaussian by default
    psf=gaussian_psf(fwhm(1)/sqrt(8*log(2)));
end

for order=orders(:).'
   if strcmp(settings.dec.psfmodel,'airy')
       PSF = map2D(airyPSF(fwhm(1)*sqrt(order)));
   else % gaussian by default
       PSF=gaussian_psf(fwhm(1)*sqrt(order)/sqrt(8*log(2)));
       PSF=PSF(:)*PSF;
   end
   [~,~,k]=size(sofi{order});
   for k=1:k
        img=abs(sofi{order}(:,:,k));
        scale=max(img(:));
        img=(scale*deconvlucy(img./scale,PSF,iter)).^lincoeffs(order);
        
        sofi{order}(:,:,k)=img;
        try  
            if settings.dec.reconvolve == 1
                sofi{order}(:,:,k)=conv2(psf,psf,img,'same');    
            end
        catch
        end
   end
end
% eof