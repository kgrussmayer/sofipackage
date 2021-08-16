function out=hriSegmentation(Id,bgth,LoGsize,out)
% Description
% ===========
% hriSegmentation filters an image with a LaplacianOfGaussian and segments
% it with a threshold for the background
%
% Input
% =====
% Id:   image to be segmented
% bgth: background threshold
%
% Output
% ======
% out.segx,out.segy:    The segments weighted center of gravity
% out.segA,out.segE:    The segments area and energy
%
% Author
% ======
%      Stefan Geissbuehler
%      Swiss Federal Institute of Technology, CH-1015 Lausanne
%      Laboratoire d'optique biomedicale, LOB
%      Biomedical imaging group, BIG

% Compute Laplacian of Gaussian
Is=LaplacianOfGaussian(Id,LoGsize);

% subtract background
out.bgmap=double((Is<bgth));
Is=out.bgmap.*Id;

% segment processed image
if any(Is(:))
    [S,n]=segmentImage(Is);

% analyze segments
    [out.segx,out.segy,out.segA,out.segE]=analyzeSegments(S,Id,n,1); % ms=1 (Minimum subtraction on)
else
    out.segx = [];
    out.segy = [];
    out.segA = [];
    out.segE = [];
end
% eof