function [imbw] = getbw(wf)

% cb = 10; % crop borders
% wf = wf(cb:end-cb+1,cb:end-cb+1); 
Ipeak = max(wf(:));        
im = uint8(wf./Ipeak*255);      

% thresholding
[imbw,~] = imthresh(im);

strElem = [1 1 1; 1 1 1; 1 1 1];

% for ii = 1:2
imbw = imerode(imbw,strElem);
% end