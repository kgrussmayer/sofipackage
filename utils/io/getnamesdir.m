function fnames = getnamesdir(folderpath, str, return_suffix)
% Find list of specific files in folder
%
% Input/output arguments:
%
%   str   ... [string]  file extension
%   fnames   ... cell   list of available files

if nargin < 3; return_suffix = false; end;

fileglob = [folderpath, filesep, '*', str];
list = dir(fileglob);
if ~isempty(list)
    num = length(list);
    for I = 1:num
        if return_suffix == true
            fnames{I} = list(I).name(1:end);
        else
            fnames{I} = list(I).name(1:end-(size(str, 2)));
        end
    end
else
	disp('No files found in specified folder. Looking in subfolders')
    list = dir([folderpath, filesep, '*',filesep,'*', str]);
    num = length(list);
    for I = 1:num
        folder=list(I).folder;
         if return_suffix == true
            fnames{I} = [folder, filesep, list(I).name(1:end)];
         else
            fnames{I} = [folder, filesep, list(I).name(1:end-(size(str, 2)))];
            
         end
    end
end
% eof