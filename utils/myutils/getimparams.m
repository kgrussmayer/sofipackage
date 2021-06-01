function [params] = getimparams(wf)
   
cb = 5; % crop borders
wf = wf(cb:end-cb+1,cb:end-cb+1); 

% wf = double(wf);
sigma1 = estimate_noise(wf);
% sigma2 = sqrt(PCANoiseLevelEstimator(wf));
    
Ipeak = max(wf(:));        
im = uint8(wf./Ipeak*255);
% sigma1 = uint8(sigma1./Ipeak*255); 

imSize = size(im);       

% thresholding
[imbw,thresh] = imthresh(im);
% thresh = graythresh(im);
% thresh = 0.5;
% [imbw] = im2bw(im, thresh);

SignalMask = imbw;
SNRbackgMask=imbw;

% create SNR background mask "far away from the signal" of the sample
strElem = [1 1 1; 1 1 1; 1 1 1];
for kk = 1:10;
SNRbackgMask = imdilate(SNRbackgMask,strElem);
% SignalMask = imerode(SignalMask,strElem);
end

SNRbackgMask = imcomplement(SNRbackgMask);
SBRbackgMask = logical(ones(imSize) - (SignalMask + SNRbackgMask));

params.SNR = 10*log10(mean(mean(wf(SignalMask).^2))/mean(mean(wf(SNRbackgMask).^2)));
% params.SNR1 = 10*log10(mean(mean(wf(SignalMask).^2))/sigma1);
% params.SNR2 = 10*log10(mean(mean(wf(SignalMask).^2))/sigma2);
params.SNR1 = 10*log10(mean(mean((wf-sigma1).^2))/sigma1);
% params.SNR2 = 10*log10(mean(mean((wf-sigma2).^2))/sigma2);

params.SBR = 10*log10(mean(mean(wf(SignalMask).^2))/mean(mean(wf(SBRbackgMask).^2)));
params.bcg = mean(mean(wf(SBRbackgMask)));
params.SB = median(median(wf(SignalMask)))/params.bcg;
params.pSB =  max(max(wf(SignalMask)))/params.bcg; %peak signal to noise rao

params.thresh = thresh;
params.sigma1 = sigma1;
% params.sigma2 = sigma2;
