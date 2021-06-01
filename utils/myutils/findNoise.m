function [noise,rect] = findNoise(im,c,display)

if nargin < 3
    display = 0;
end
if nargin < 2
    c = 50;
end

if c < 10
    c = 10;
end

nPix = length(im);

newPix = c.*floor(nPix/c);
im = im(1:newPix,1:newPix);

% divide the image into block of 20 pix
% and compute the local contrast Imax-Imin

cMap = [];
cx = 1;
cy = 1;
for k = 1:c:newPix
    for j = 1:c:newPix
        temp = im(k:k+c-1,j:j+c-1);
        cMap(cy,cx) = (max(temp(:))-min(temp(:)));
        cx = cx + 1;
    end
    cy = cy + 1;
    cx = 1;
end

% pic the location with lowest contrast
[~,ind] = min(cMap(:));
[y,x] = ind2sub(size(cMap),ind);

indy = ((y-1)*c+1):(y*c)-5;
indx = ((x-1)*c+1):(x*c)-5;

noise = im(indy,indx);
rect = [indx(1) indy(1) c-5 c-5].*nPix/newPix;

if display
    figure(display)
    subplot(221);imagesc(im)
    subplot(222);imagesc(cMap)
    subplot(223);imagesc(noise)
    subplot(224);imagesc(im);hold on
    rectangle('position',[indx(1) indy(1) c-5 c-5],'edgecolor','w')
    hold off
end
