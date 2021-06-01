function [sofi]=sofiLinearize4(sofi,fwhm,noise,orders,iter,lincoeffs,settings)
%Deconvolve and denoise flat cumulants and linearize the brightness
%by taking the lincoeffs-th root of the deconvolved cumulants. The result
%is reconvolved with a realistic point-spread function.
%
%Inputs:
% sofi          Flat cumulants (see sofiFlatten)
% fwhm          PSF diameter (see sofiFlatten)
% noise         Noise level                         {1%}
% orders        Cumulant orders                     {all}
% iter          Maximum deconvolution iterations    {10}
% lincoeffs     Adaptive linearization coefficients
%
%Output:
% sofi          Linarized cumulant images

if nargin < 3 || isempty(noise)
   noise=0.01;
end

if nargin < 4 || isempty(orders)
   orders=find(~cellfun(@isempty,sofi));
end

if isempty(orders)
   return;
end
if nargin < 5 || isempty(iter)
   iter=10;
end
% imgpar2=cell(1,length(orders));
% imgpar3=cell(1,length(orders));

%fwhm=fwhm(1)/sqrt(8*log(2));
%psf=gaussian(fwhm);

if strcmp(settings.dec.psfmodel,'airy')
    psf = airyPSF(fwhm(1));
else % gaussian by default
    psf=gaussian(fwhm(1)/sqrt(8*log(2)));
end


for order=orders(:).'
   if strcmp(settings.dec.psfmodel,'airy')
       PSF = map2D(airyPSF(fwhm(1)*sqrt(order)));
   else % gaussian by default
       PSF=gaussian(fwhm(1)*sqrt(order)/sqrt(8*log(2)));
       PSF=PSF(:)*PSF;
   %PSF=gaussian(fwhm*sqrt(order));
   %PSF=PSF(:)*PSF;
   end
   [~,~,k]=size(sofi{order});
   for k=1:k
        img=abs(sofi{order}(:,:,k));
        scale=max(img(:));
        img2=(scale*deconvlucy(img./scale,PSF,iter));
      
%       gamma= gammaest2(img,img2,'mean');
%       gamcor(order) = lincoeffs(order)/gamma;
        img3=(img2.^lincoeffs(order));
        
        sofi{order}(:,:,k)=img3;
        try  
            if settings.dec.reconvolve == 1
                sofi{order}(:,:,k)=conv2(psf,psf,img3,'same');    
            end
        catch
            
        end      
        

   end
end
