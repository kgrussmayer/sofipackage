function [X_shift, Y_shift] = ccrShiftEstimation(im1,im2,upsc)
% Estimate shift in between two images (img1 and img2) in x and y direction
% using cross corelation upsc is the upscaling factor, it upscales the
% precision to obtain subpixel shifts. Final precision is
% (- 1/upsc, + - 1/upsc). Due to computantional demands, upsc bigger then
% 20 can lead to lack of the memory
%
%Inputs:
% im1       first image
% im2       second image
% upsc      upscaling factor
%
%Outputs:
% X_shift   shift in x direction (image columns)
% Y_shift   shift in y direction (image rows)

% Copyright © 2018 Tomas Lukes
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

if size(im1,3) > 1 % check if the images are color images
    im1 = rgb2ycbcr(im1);
    img1 = double(im1(:,:,1));
else
    img1 = double(im1);
end

if size(im2,3) > 1
    im2 = rgb2ycbcr(im2);
    img2 = double(im2(:,:,1));
else
    img2 = double(im2);
end

img1FT = fft2(img1);
img2FT = fft2(img2);

% upscaling factor   - > upscale the precision,
[rows,cols]=size(img1FT);
mlarge=rows*upsc;
nlarge=cols*upsc;

% zero padding
CC=zeros(mlarge,nlarge);

% shift the fft to the center
CC(rows+1-fix(rows/2) : rows+1+fix((rows-1)/2), ...
    cols+1-fix(cols/2) : cols+1+fix((cols-1)/2)) = ...
    fftshift(img1FT).*conj(fftshift(img2FT));

% compute crosscorrelation
CC = ifft2(ifftshift(CC));
% find the peak
peakMax = max(CC(:));
[peakY,peakX] = find(CC == peakMax);

% get the shift in the original pixel grid from the position of the
% crosscorrelation peak
[rows,cols] = size(CC);

if peakY > fix(rows/2)
    Y_shift = peakY - rows - 1;
else
    Y_shift = peakY - 1;
end
if peakX > fix(cols/2)
    X_shift = peakX - cols - 1;
else
    X_shift = peakX - 1;
end

% calculate the output shift
Y_shift=Y_shift/upsc;
X_shift=X_shift/upsc;
% eof