% Compute the normalized cross-correlation of the signal 1 and 2
% cc = normxcorr(sig1,sig2)

function [cc, lags] = normxcorr(refSig,measSig)

l1 = length(refSig);
l2 = length(measSig);

% suppress the mean
refSig = refSig-mean(refSig);
measSig = measSig-mean(measSig);
% cross-correlate
if nargout>1
    [cc, lags] = xcorr(measSig,refSig);
else
    cc = xcorr(measSig,refSig);
end
% normalize
cc = sqrt(l2*l1).*cc./(sqrt(sum(refSig.*refSig)*sum(measSig.*measSig)).*l2);
