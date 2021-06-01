function show2subs(im1,im2,figID)
if nargin < 3; figure;
else
    figure(figID)
end
subplot(121);imshow(im1,[])
subplot(122);imshow(im2,[])
