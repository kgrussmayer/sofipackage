%sofi=sofiLinearize(sofi,fwhm,noise,orders,iter)
%-----------------------------------------------
%
%Deconvolve and denoise flat cumulants and linearize the brightness
%by taking the order-th root of the deconvolved cumulants. The result
%is reconvolved with a realistic point-spread function.
%
%Inputs:
% sofi      Flat cumulants (see sofiFlatten)
% fwhm      PSF diameters (see sofiFlatten)
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
function sofi=sofiLinearize3D_TL(sofi,fwhm,noise,orders,iter)
if nargin < 3 || isempty(noise)
   noise=0.01;
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
dims=numel(fwhm);
fwhm=fwhm/sqrt(8*log(2));
psf=gaussian(fwhm(1));
psf=psf(:)*gaussian(fwhm(2));
if dims == 3
   [x,y]=size(psf);
   psf=psf(:)*gaussian(fwhm(3));
   psf=reshape(psf,x,y,size(psf,2));
end
for order=orders(:).'
   sigma=fwhm*sqrt(order);
   PSF=gaussian(sigma(1));
   PSF=PSF(:)*gaussian(sigma(2));
   images=abs(sofi{order});
   sizes=size(images);
   if dims == 3
      [x,y]=size(PSF);
      PSF=PSF(:)*gaussian(sigma(3));
      PSF=reshape(PSF,x,y,size(PSF,2));
      k=prod(sizes(4:end));
   else
      k=prod(sizes(3:end));
      images=reshape(images,[sizes(1:2) 1 k]);
   end
   for k=1:k
      image=images(:,:,:,k);
      scale=max(image(:));
      image=(scale*deconvlucy(image/scale,PSF,iter)).^(1/order);
      image(image < order*noise*max(image(:)))=0;
      images(:,:,:,k)=convn(image,psf,'same');
   end
   sofi{order}=reshape(images,sizes);
end
