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

%Copyright � 2012 Marcel Leutenegger et al, �cole Polytechnique F�d�rale de Lausanne,
%Laboratoire d'Optique Biom�dicale, BM 5.142, Station 17, 1015 Lausanne, Switzerland.
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
function sofi=sofiLinearize(sofi,fwhm,noise,orders,iter,settings)
if nargin < 3 || isempty(noise)
   noise=0.01;
%    noise = 0.1;
end
if nargin < 4 || isempty(orders)
   orders=find(~cellfun(@isempty,sofi));
end
orders=orders(orders > 1);
if isempty(orders)
   return;
end
if nargin < 5 || isempty(iter)
   iter=10;
end
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
    end

   [~,~,k]=size(sofi{order});
   for k=1:k
      img=abs(sofi{order}(:,:,k));
      scale=max(img(:));
      img=(scale*deconvlucy(img/scale,PSF,iter)).^(1/order);
      
      sofi{order}(:,:,k)=img;
       
      try
          if settings.dec.reconvolve == 1
            sofi{order}(:,:,k)=conv2(psf,psf,img,'same');    
          end
      catch
      end

   end
end
