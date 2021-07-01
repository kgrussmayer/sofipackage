function unames = uniNames(fnames,fnames2)

names = unique(fnames2);

for ii = 1:numel(names)
    count = 1;
    for jj = 1:numel(fnames2)
        if strcmp(names{ii},fnames2{jj})
            unames(ii).fnames{count} = fnames{jj};
            unames(ii).fnames2{count} = fnames2{jj};
            unames(ii).fnum{count} = jj;
            count = count +1;
        end
        
    end
end
% eof