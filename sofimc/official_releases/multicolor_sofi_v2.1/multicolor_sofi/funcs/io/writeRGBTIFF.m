function writeRGBTIFF(data, tiffname,res)
% Write uint8 RGB TIF stack
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

filename = [tiffname '.tiff'];
maxCount = 10;
count = 1;
while exist(filename, 'file') == 2
    filename = [tiffname '_' num2str(count) '.tiff'];
    count = count +1;
    if count == maxCount
        break;
    end
end

    for k = 1:size(data,4)
        t = Tiff(filename, 'a');
        tagstruct.ImageLength = size(data, 1);
        tagstruct.ImageWidth = size(data, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.XResolution = res(1);
        tagstruct.YResolution = res(2); % YResolution encode z-res
        tagstruct.Photometric = 1;
        tagstruct.BitsPerSample =  8;                        % float data
        tagstruct.SamplesPerPixel = 3;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        t.setTag(tagstruct);
        t.write(data(:,:,:,k));
        t.close();
    end

end