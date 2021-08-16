function fnames = getnames(folderpath, str)
% Find list of specific files in folder
%
% Input/output arguments:
%
%   str   ... [string]  mask of functions
%   fnames   ... cell   list of available files

list = dir([folderpath,filesep,str]);
num = length(list);

for I = 1:num
    fnames{I} = list(I).name(1:end-4);
end

if isempty(I)
    fnames = [];
end
% eof