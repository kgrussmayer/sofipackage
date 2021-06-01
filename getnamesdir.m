function fnames = getnamesdir(folderpath, str)
% Find list of specific files in folder
%
% Input/output arguments:
%
%   str   ... [string]  file extension
%   fnames   ... cell   list of available files

% fd = fileparts_dir([mfilename('fullpath') '.m']);
% folderpath = 'H:\tlukes\main_simulation_v05\outputsim\sim1';
% str = '';
fileglob = [folderpath, filesep, '*', str]
list = dir(fileglob)
if ~isempty(list)
    num = length(list);
    for I = 1:num
      fnames{I} = list(I).name(1:end-(size(str, 2)));
    %   fnames{I} = list(I).name(1:end);
    end
else
	disp('No files found in specified folder. Looking in subfolders')
    list = dir([folderpath,filesep, '*',filesep,'*', str]);
    num = length(list);
    for I = 1:num
        folder=list(I).folder(size(folderpath,2)+2:end);
        fnames{I} = [folder, filesep, list(I).name(1:end-(size(str, 2)))];
    %   fnames{I} = list(I).name(1:end);
    end
end