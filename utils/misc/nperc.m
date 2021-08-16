function im = nperc(im,np)
% saturate np % of highest image values

% find values within the percentage (highest n % values)
vs=sort(im(:),'descend');
n=round(numel(im(:))*np/100)+1;

% the result
thresh=min(vs(1:n));

im = min(im,thresh);
% eof