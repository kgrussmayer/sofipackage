function [params] = getimparams(wf)

cb = 5; % crop borders
wf = wf(cb:end-cb+1,cb:end-cb+1);

sigma1 = estimate_noise(wf);

Ipeak = max(wf(:));
im = uint8(wf./Ipeak*255);

imSize = size(im);

% thresholding
[imbw,thresh] = imthresh(im);

SignalMask = imbw;
SNRbackgMask=imbw;

% create SNR background mask "far away from the signal" of the sample
strElem = [1 1 1; 1 1 1; 1 1 1];
for kk = 1:10
    SNRbackgMask = imdilate(SNRbackgMask,strElem);
end

SNRbackgMask = imcomplement(SNRbackgMask);
SBRbackgMask = logical(ones(imSize) - (SignalMask + SNRbackgMask));

params.SNR = 10*log10(mean(mean(wf(SignalMask).^2))/mean(mean(wf(SNRbackgMask).^2)));
params.SNR1 = 10*log10(mean(mean((wf-sigma1).^2))/sigma1);

params.SBR = 10*log10(mean(mean(wf(SignalMask).^2))/mean(mean(wf(SBRbackgMask).^2)));
params.bcg = mean(mean(wf(SBRbackgMask)));
params.SB = median(median(wf(SignalMask)))/params.bcg;
params.pSB =  max(max(wf(SignalMask)))/params.bcg; %peak signal to noise rao

params.thresh = thresh;
params.sigma1 = sigma1;
% eof