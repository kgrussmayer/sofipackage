function im_rgb = mergeToRgb(im_red, im_green, im_blue)
% Create RGB image.
%
% Inputs:
% im_red     image corredponding to the red channel (rows, cols)
% im_green   image corredponding to the green channel (rows, cols)
% im_blue    image corredponding to the blue channel (rows, cols)
% 
% Outputs:
% im_rgb    RGB image (rows, cols, RGB channels)

% Copyright © 2018 Tomas Lukes
% École Polytechnique Fédérale de Lausanne,
% Laboratory of Nanoscale Biology, http://lben.epfl.ch/
 
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

im_blue=mean(im_blue,3); 
im_blue(im_blue<0) = 0;
im_blue_norm=im_blue/max(im_blue(:));

im_green=mean(im_green,3); 
im_green(im_green<0) = 0; 
im_green_norm=im_green/max(im_green(:));

im_red=mean(im_red,3); 
im_red(im_red<0) = 0;
im_red_norm=im_red/max(im_red(:));

im_rgb =cat(3,im_red_norm,im_green_norm,im_blue_norm);
% eof