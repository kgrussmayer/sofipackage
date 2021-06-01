function saveSettingsTxt(settings,fname)
% Save settings to txt file.
%
% Inputs:
% settings      [struct] settings used for the processing
% fname         name of the file to be saved

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

fList = fieldnames(settings);
fid = fopen( fname, 'wt' );

for k = 1:length(fList)
    
    if ~isstruct(getfield(settings,fList{k}))
        t = getfield(settings,fList{k});
        if isstr(t)
            fprintf(fid,'%s',['settings.',fList{k},' = ',t]);
        else
            fprintf(fid,'%s',['settings.',fList{k},' = ',num2str(t)]);
        end
        
            fprintf(fid,'\n');
    else
        fprintf(fid,'\n');
        subStruct = getfield(settings,fList{k});
        subfList = fieldnames(subStruct);
            for n = 1:length(subfList)
                t = getfield(subStruct,subfList{n});
                	if isstr(t)
                        fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ',t]);
                    else
                        fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ',num2str(t)]);
                    end
                    fprintf(fid,'\n');
            end
    end
    
end

fclose(fid);