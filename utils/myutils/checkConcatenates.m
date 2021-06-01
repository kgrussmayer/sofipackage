function stack
% find numbers, indentify groups of files with the same names
fnames2 = fnames;
clear t;

for ii = 1:numel(fnames)
    str = [fnames{ii},''];
    strnum = regexp(str,'(\d*).ome','tokens','once');
    
    if ~isnan(str2double(strnum))
        t(ii) = str2double(strnum);
        str2 = str(1:end-length(strnum{1})-4);
        fnames2{ii} = str2;
        
    else
        t(ii) = NaN;
        ii
    end

end

% unames2 = uniNames(fnames,fnames2);

%%
fnums = 1:length(t);
fnums(isnan(t)) = []; % stores order of the fif files in fnames array
t(isnan(t)) = []; %stores number of the tif files

[t, indx] = sort(t);
fnums = fnums(indx);

p=find(diff(t)==1); % TODO najit v diff(t) segmenty po sobe jdoucich jednice, u tehle kontrola jestli maji stejne jmeno, pokud ano, tak sloucit
q=[p,p+1];
q = unique(q);
p = fnums(q);% positions of all the pairs of consecutive numbers
 

unames = uniNames(fnames(p),fnames2(p));


ngroups = numel(unames);
stack = [];
settings.io.concTifSave =1;

for ii = 1:ngroups
    numbers = cell2mat(unames(ii).fnum);
    [numbers, order] = sort(numbers);
    
    names = unames(ii).fnames;
    names = names(order);
    
    % save big tif file
    if settings.io.concTifSave == 1 
        
    % load tif files
    for jj = 1:numel(names)
        settings.io.imageName = names{jj};
        settings.io.imageFile = [settings.io.imagePath, filesep,settings.io.imageName];
        [temp,frames] = loadStack(settings);
        stack = cat(3,stack,temp);
    end
   
%         options.big = true; % Use BigTIFF format
%         options.compress = 'no';
%         tic
%         saveastiff(stack, 'BigTiff(2GB+2GB).btf', options);
%         toc
        
%         tic
%         save('testBig.mat','-v7.3','stack');
%         toc

        tic; 
        saveFastMat('testBig_savefast2',stack);
        toc
        
    end

end