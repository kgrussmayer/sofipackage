function im_rgb = mergeToRgb(im_red, im_green, im_blue)

% normalize all channels with the same value
nCoef = max([max(im_blue(:)) max(im_green(:)) max(im_red(:))]);

im_blue=mean(im_blue,3); 
im_blue(im_blue<0) = 0;
im_blue_norm=im_blue/max(im_blue(:));

im_green=mean(im_green,3); 
im_green(im_green<0) = 0; 
im_green_norm=im_green/max(im_green(:));

im_red=mean(im_red,3); 
im_red(im_red<0) = 0;
im_red_norm=im_red/max(im_red(:));

im_rgb =cat(3,im_red_norm,im_green_norm,im_blue_norm);