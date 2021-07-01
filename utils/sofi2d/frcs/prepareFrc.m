function [im1, im2] = prepareFrc(sofi_lin,settings,order)
% prepare data for FRC calculation

im1 = mean(sofi_lin{order}(:,:,1:2:end),3);
im2 = mean(sofi_lin{order}(:,:,2:2:end),3);
% eof
