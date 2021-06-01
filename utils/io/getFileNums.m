function [fnames2,t] = getFileNums(fnames)
% detect and cut numbers at the end of input image names

fnames2 = fnames;

for ii = 1:numel(fnames)
    str = [fnames{ii},''];
    strnum = regexp(str,'(\d*).ome','tokens','once');
    
    if ~isnan(str2double(strnum))
        t(ii) = str2double(strnum);
        str2 = str(1:end-length(strnum{1})-4);
        fnames2{ii} = str2;
        
    else
        t(ii) = NaN;
    end

end