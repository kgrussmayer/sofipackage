function [stack, settings] = concatFiles(fnames,settings, str)
% find numbers, indentify groups of files with the same names

[fnames2,t] = getFileNums(fnames);

fnums = 1:length(t);
fnums(isnan(t)) = []; % stores order of the fif files in fnames array
t(isnan(t)) = []; %stores number of the tif files

[t, indx] = sort(t);
fnums = fnums(indx);

p=find(diff(t)==1);
q=[p,p+1];
q = unique(q);
p = fnums(q);% positions of all the pairs of consecutive numbers

if isempty(p)
    p=1;
    q=1;
end

unames = uniNames(fnames(p),fnames2(p));

stack = [];

numbers = cell2mat(unames(1).fnum);
[numbers, order] = sort(numbers);

names = unames(1).fnames;
names = names(order);

% load tif files
for jj = 1:numel(names)
    settings.io.imageName = names{jj};
    settings.io.imageFile = [settings.io.imagePath,filesep,settings.io.imageName];
    [temp,frames] = loadStack(settings, str);
    stack = cat(3,stack,temp);
end

settings.io.imageName = unames(1).fnames2{1};
settings.io.imageFile = [settings.io.imagePath,filesep,settings.io.imageName];

[~, settings.io.imageName,~]=fileparts(settings.io.imageName);

% save big tif file
if settings.io.concatSave == 1 && exist([settings.io.imageFile,'.mat'], 'file') ~= 2
    
    tic;
    saveFastMat(settings.io.imageFile,stack);
    toc
    
end
% eof