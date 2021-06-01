function imrgb = makeComposite(im1,im2)
im1 = double(im1); im2 = double(im2);
imrgb(:,:,1) = im1./max(im1(:));
imrgb(:,:,2) = im2./max(im2(:));
imrgb(:,:,3) = im2./max(im2(:));