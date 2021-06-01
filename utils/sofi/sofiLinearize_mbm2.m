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

function [sofi2]=sofiLinearize_mbm2(sofi,fwhm,orders,iter,lincoeffs,fcpx)
   

for order=orders(:).'
    fwhm1=(fwhm)/sqrt(8*log(2));
%     PSF=gaussian(fwhm1/sqrt(order));
    PSF=gaussian(fwhm1*sqrt(order));
%     PSF=gaussian(fwhm1);
    PSF=PSF(:)*PSF;
    [~,~,k]=size(sofi{order});

    for k=1:k
        img=sofi{order}(:,:,k);

        [img] = mseBregDeconv3(img,PSF,iter,fcpx,[],[],lincoeffs(order));
        
%         [img] = mseBregDeconv3b(img,PSF,iter,fcpx,[],[],[]);
%         img=sofiBregman2D(img,fwhm,order,iter,[],0.01,1e-6);
        
        sofi2{order}(:,:,k)=img;
   end
end
