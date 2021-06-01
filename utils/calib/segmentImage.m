%[S,n]=segmentImage(I)
%---------------------
%
%Segment an image I and return the segment map S and the number of segments n.
%
function [I,n]=segmentImage(I)
I(end+1,end+1)=0;
%
% Segment along x
%
a=find(I);
I(a)=cumsum([1;diff(a) > 1]);
I=I.';
%
% Segment along y
%
a=find(I);
b=I(a);
for n=find(diff(a) == 1).'
   b(b == b(n+1))=b(n);
end
%
% Enumerate segments
%
c=unique(b);
n=numel(c);
d=zeros(c(end),1);
d(c+0)=1:n;          % avoid logicals
I(a)=d(b+0);
I=I(1:end-1,1:end-1).';
