function stack = load_tifFile(file_name, n_frames,roi)
% Load a stack of images from a multiple tif file.
% Fname contains path to the file to be loaded. Nframes is a number of
% frames which are going to be loaded from the multiple tif file. Optinally
% only roi (region of interest) can be loaded.
%
% Inputs:
% file_name     path to the multiple tif file
% n_frames      number of frames to be loaded
% roi           region of interest to be loaded
%
% Outputs:
% stack         stack of images (rows, columns, frames)

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

fileinfo = imfinfo(file_name);

if nargin < 3
    roi = [];
end

if nargin < 2 || ~any(n_frames) || isempty(n_frames)
    n_frames=length(fileinfo);
end

if n_frames > length(fileinfo)
    n_frames=length(fileinfo);
    disp(['Tiff file contains only ',num2str(n_frames),' images'])
end

fig=statusbar('Loading data...');

for ii=1:n_frames
    fig=statusbar(ii/n_frames,fig);
    temp=uint16(imread(file_name,ii,'info',fileinfo));
    if size(temp,3) > 1 % check RGB image
        temp = rgb2gray(temp);
    end
    
    if any(roi)
        stack(:,:,ii) = temp(roi(1,:),roi(2,:));
    else
        stack(:,:,ii) = temp;
    end
end
delete(fig);
% eof