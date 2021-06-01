function [corela, corelf] = stackcorel(stack,MaxCorrSamp)
% CALCULATE AUTOCORRELATION ON IMAGE STACK
% IN
% stack ... Image Stack.
% samples ... number of samples of the autocorelation funciton
%
% OUT
% corela = xy average autocorrelation over time

stack = stack(:,:,1:MaxCorrSamp);
[sy,sx,sz] = size(stack);

mask = var(single(stack),[],3);
mask = mean(single(stack),3);
mask = mask./max(mask(:));
% mask(mask>0.8) = 0;
% mask = mask./max(mask(:));
mask_f = logical(repmat(im2bw(mask,graythresh(mask)),1,1,sz));
mask_b = logical(repmat(imcomplement(im2bw(mask,0.8*graythresh(mask))),1,1,sz));

corela = abs(ifft(abs(fft(cat(3,stack,zeros(sy,sx,sz)),[],3)).^2,[],3));
corela = corela(:,:,1:sz);

corelf = corela;
corelf(mask_f) = [];
corelf = reshape(corelf,[],sz);
corelf = max(corelf,[],1);
corelf = corelf(1:MaxCorrSamp);
corelf = corelf./max(corelf(:));

% corelb = corela;
% corelb(mask_b) = [];
% corelb = reshape(corelb,sz,[]);
% corelb = mean(corelb,2);
% corelb = corelb(1:MaxCorrSamp);
% corelb = corelb./max(corelb(:));

corela = squeeze(mean(mean(corela,1),2));
corela = corela(1:MaxCorrSamp);
corela = corela./max(corela(:));
end

% %%
% figure,
% plot(corela);
% hold on;
% plot(corelb);
% plot(corelf);