function saveFigure(fh,output_path,file_name,format,add_id)
% Save the figure fh according to output_path and file_name
%
% Inputs:
% fh                figure handle
% output_path       output path
% file_name         name of the output file
% format            format of the output file i.e. pdf, png, tif, fig
% add_id             add a unique identifier to the file name yes/no {1, 0}
%
% Example of use:
% save the current figure in tif file
% fh = gcf;saveFigure(fh,outputPath,'nameOfTheFigure','tif');

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
if nargin <4
    format='fig';
    add_id = 0;
end

if nargin <5
    add_id = 0;
end

% create the output directory if it does not exists
if (~exist(output_path,'dir'))
    mkdir(output_path);
end

if add_id == 1
    nowID = datestr(now,'mmmm-dd-yyyy_HHMMSS');
    % add a unique identifier to the file name
    outputFile = [output_path,filesep,file_name,'_',nowID];
else
    outputFile = [output_path,filesep,file_name];
end

set(fh, 'PaperPositionMode','auto');
saveas(fh,[outputFile,'.',format]);
disp(['figure was succesfully saved into: ',outputFile]);
% eof