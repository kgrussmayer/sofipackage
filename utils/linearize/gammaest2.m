function gamma = gammaest2(img1,img2,type)
% estimate gamma between two images
    
% img1 = img1./max(img1(:));
% 
% img2 = img2./max(img2(:));

% q = min(img1(img1>0))*1e-40;
% q=+1e-16; % min value to prevent zero values
img1 = imfilter(abs(img1),fspecial('gaussian',[9 9],1));
img2 = imfilter(abs(img2),fspecial('gaussian',[9 9],1));

q = 0;
img1 = img1(:);
img2 = img2(:);
mask = (img1>0 & img2>0);
img1(mask==0)=[];
img2(mask==0)=[];

denom = (log10(img1(:)+q));
gamma = (log10(img2(:)+q))./denom;
gamma(denom==0)=[];

if strcmp(type,'mean')
    gamma = mean(gamma);
elseif strcmp(type,'median')
    gamma = median(gamma);
else
    disp('type of calculation was not recognized');    
end

