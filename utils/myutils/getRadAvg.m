function r = getRadAvg(im)

if size(im,1) ~= size(im,2)
    error('getRadAvg supports only square input');
elseif length(size(im)) ~= 2
    error('getRadAvg supports only 2D matrix as input');
end

r = mean(im2pol(im),1);
r = imresize(r,[1 ceil(size(im,2)/2)],'bilinear');