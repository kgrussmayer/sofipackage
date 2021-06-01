function [ XX ] = stackcorel(stack,MaxCorrSamp)
% CALCULATE AUTOCORRELATION ON IMAGE STACK
% IN
% stack ... Image Stack.
% samples ... number of samples of the autocorelation funciton
%
% OUT
% corela = xy average autocorrelation over time

[sy,sx,sz] = size(stack);

corela = abs(ifft(abs(fft(cat(3,stack,zeros(sy,sx,sz)),[],3)).^2,[],3));
corela = squeeze(mean(mean(corela,1),2));
XX.Corela = corela(1:MaxCorrSamp);

end