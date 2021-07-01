function [imbw] = getbw(wf)

Ipeak = max(wf(:));        
im = uint8(wf./Ipeak*255);      

% thresholding
[imbw,~] = imthresh(im);

strElem = [1 1 1; 1 1 1; 1 1 1];

imbw = imerode(imbw,strElem);
% eof