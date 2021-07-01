function writeTIFF(data,tiffname,res)
% Write image stack as multi frame 32bit TIF
%
% Inputs:
% data          path to the multiple tif file
% tiffname      name of the file to be saved
% res           image resolution (optional)

% Copyright © 2018 Adrien Descloux
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


if nargin < 3
    res = [1 1]; % save generic resolution if none provided
end

filename = [tiffname '.tif'];
maxCount = 10;
count = 1;
while exist(filename, 'file') == 2
    filename = [tiffname '_' num2str(count) '.tif'];
    count = count +1;
    if count == maxCount
        break;
    end
end

for k = 1:size(data,3)
    t = Tiff(filename, 'a');
    tagstruct.ImageLength = size(data, 1);
    tagstruct.ImageWidth = size(data, 2);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.XResolution = res(1);
    tagstruct.YResolution = res(2); % YResolution encode z-res
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = 1;
    tagstruct.BitsPerSample =  32;                        % float data
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tagstruct);
    t.write(single(data(:,:,k)));
    t.close();
    pause(0.01)
end
% eof