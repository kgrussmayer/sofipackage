function imgs = interpolateDriftcor(imgs, drifts)

[rows, columns, frames] = size(imgs);
[xx, yy]=meshgrid(1:columns,1:rows);

for k=1:frames
    disp(['dcor: ',num2str(k)]);
    img = imgs(:,:,k);

    imReg=interp2(xx,yy,img,xx-drifts(k,1),yy-drifts(k,2),'linear');
    imReg(isnan(imReg)) = 0;
    imgs(:,:,k)=imReg;
end
% eof