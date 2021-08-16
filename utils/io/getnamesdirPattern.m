function fnames = getnamesdirPattern(folderpath, pattern)
% Find list of specific files in folder
% Tomas Lukes, tomas.lukes@epfl.ch
%
% Input/output arguments:
%
%   folderpath  .. path with folders of the measured data
%   fnames   ... cell   list of available files

if ~isempty(pattern)
    list = dir(fullfile(folderpath, pattern));
else
    list = dir(folderpath);
end
isub = [list(:).isdir]; 

fnames = {list(isub).name}';
fnames(ismember(fnames,{'.','..','calibration'})) = [];

disp('Data to be processed:');
disp(fnames);
% eof