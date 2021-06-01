function plotImagesAll(imgs,names,cmap,brrange,roi)



if nargin < 4
    brrange = [0 1];
end

N = size(imgs,3);


for n = 1: N

subplot_tight(ceil(N/ceil(sqrt(N))),ceil(sqrt(N)), n, [.02 .01])

imshow(imgs(:,:,n),brrange);colormap(cmap);hold on;
title(['im',num2str(n),' : ',names{n}]);colormap(cmap);

if exist('roi','var')
    plotRect(roi,'R1','red');
end

end