%sofi=sofiLinearize(sofi,fwhm,noise,orders,iter)
%-----------------------------------------------
%
%Deconvolve and denoise flat cumulants and linearize the brightness
%by taking the order-th root of the deconvolved cumulants. The result
%is reconvolved with a realistic point-spread function.
%
%Inputs:
% sofi      Flat cumulants (see sofiFlatten)
% fwhm      PSF diameter (see sofiFlatten)
% noise     Noise level                         {1%}
% orders    Cumulant orders                     {all}
% iter      Maximum deconvolution iterations    {10}
%
%Output:
% sofi      Linar cumulant images

%Copyright © 2012 Marcel Leutenegger et al, École Polytechnique Fédérale de Lausanne,
%Laboratoire d'Optique Biomédicale, BM 5.142, Station 17, 1015 Lausanne, Switzerland.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [sofi2]=sofiLinearize_mb(sofi,fwhm,orders,iter,lincoeffs,fcpx,nss)


%%% deconv settings
settings.reltol = 1e-5;
settings.maxIter = iter;
settings.numIm = nss;
settings.coeff_A = 0.0;
apodize = 1;      
wlen = nss;
% fwhm = 1.5;
for order=orders(:).'
    fwhm1=(fwhm)/sqrt(8*log(2));
    PSF=gaussian(fwhm1*sqrt(order));
    PSF=PSF(:)*PSF;
    [~,~,k]=size(sofi{order});

    for k=1:floor(k/wlen)
        img=(sofi{order}(:,:,1+(k-1)*wlen:k*wlen)); % multi deconvolution
        ming = mean(img,3);
        scale=max(img(:));
        
        lambda = estimate_noise(mean(img,3)./max(ming(:)));
        if lambda < 0.01; lambda = 0.001; end;
        if lambda > 0.1; lambda = 0.1; end;
%         uplim = 0.1;
%         if 3*sigmaNoise > uplim; settings.lambda2 = 2*sigmaNoise;
%         elseif 2*sigmaNoise > uplim; settings.lambda2 = 1*sigmaNoise;
%         elseif sigmaNoise > uplim; settings.lambda2 = uplim;
%         else settings.lambda2 = 3*sigmaNoise;
%         end
        
%         settings.lambda2 = 0.1;
        
        if order >1; omega = 0.7;else omega = 1; end
        [sy, sx] = size(img);
        Aaray = [settings.maxIter,apodize,0.04,omega,settings.numIm];
        img = img/scale;
%         img = img.^lincoeffs(order);
        img=imtophat(img,strel('disk',round(sqrt(sx^2+sy^2)/10)));
        img2 = SOFI_X7(img,PSF,Aaray);

%         gamma = gammaest2(abs(ming/scale),abs(img2),'mean');
%         gamma = 1;
        img2 = abs(img2);

        img2 = img2./max(img2(:));
        img2 = imadjust(img2,[lambda, 1],[0 1], 1);
        img2 = scale*img2;

%         img2=img2.^lincoeffs(order);
% %         img2=img2.^lincoeffs(order);

%         disp(['gamcor',num2str((1/gamma)*lincoeffs(order))])

        sofi2{order}(:,:,k)=img2;

   end
end
