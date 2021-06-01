function [im1, im2] = prepareFrc(sofi_lin,settings,order)
% prepare data for FRC calculation

im1 = mean(sofi_lin{order}(:,:,1:2:end),3);
im2 = mean(sofi_lin{order}(:,:,2:2:end),3);

% sc_im1 = max(im1(:));
% sc_im2 = max(im2(:));
% im1 = im1./sc_im1;
% im2 = im2./sc_im2;
% 
% % in1 = imadjust(in1,[thresh 1],[0 1]);
% % in2 = imadjust(in2,[thresh 1],[0 1]);
% 
% if settings.frc.bcgsub > 0 && order > 1
%     im1 = im1 - bcgsub*median(im1(:)); 
%     im2 = im2 - bcgsub*median(im2(:));
%     im1(im1<0) = 0;
%     im2(im2<0) = 0;
% end
     
% bcg1 = median(im1(:));
% bcg2 = median(im2(:));
% im1 = imadjust(im1,[settings.bcgsub*bcg1 1],[0 1]);
% im2 = imadjust(im2,[settings.bcgsub*bcg2 1],[0 1]);
