%%% run multicolor sofi
%%% lukestom 19.7. 2014
% v. 26.7.2014 - at the beginning load path to multiple tiff files
% v. 30.7.2014 - load bin files

close all; clear; 

addpath(genpath('funcs'));

%
%%% Initialization
%
addpath('funcs');
% addpath('cornerDetector');
% path to calibration files 
% file1 = 'C:\Users\Lukestom\Documents\ÈVUT\Ph.D. III.semestr\SOFI method\DATA\7\calibration\001_DU897_BV_1936.tif';
% file2 = 'C:\Users\Lukestom\Documents\ÈVUT\Ph.D. III.semestr\SOFI method\DATA\7\calibration\001_Luc247_MONO_0655.tif';


%%% load the calibration files

%
%%% TIFF FILES
%

% file1='K:\sharipov\20140529_Hela_3color\Calibration\001_DU897_BV_1936.tif';
% file2='K:\sharipov\20140529_Hela_3color\Calibration\New folder\001_Luc247_MONO_0655.tif';

% Nframes = 1;
% imc1 = load_tifFile(file1, Nframes);
% imc2 = load_tifFile(file2, Nframes);
% imc2 = flipud(imc2);
% 
% imc1(imc1==0)=min(imc1(imc1>0));
% imc2(imc2==0)=min(imc2(imc2>0));

%
%%% BIN FILES
%

file1='H:\Sharipov\Measurements\20140728_Hela\Calibration\001\data1.bin';
file2='H:\Sharipov\Measurements\20140728_Hela\Calibration\001\data2.bin';

fn1 = 'data1.bin';
pn1 = 'H:\Sharipov\Measurements\20140728_Hela\Calibration\001\';
fn2 = 'data2.bin';
pn2 = pn1;

% Nframes = 100;
imc1 = loadbinfile(pn1,fn1);
imc2 = loadbinfile(pn2,fn2);
[~,indcom]=max(squeeze(mean(max(imc1))).*squeeze(mean(max(imc2))));
im_fix=imc1(:,:,indcom);
im_mov=imc2(:,:,indcom);
% im_mov=flipdim(im_mov,2);
% im_mov=flipdim(im_mov,1);

show2subs(im_fix,im_mov)

%
%%% Find feature points
%

compTransfMask_TL(im_fix,im_mov)
% save('chosenPoints20140728_Hela.mat','pts1','pts2','tform')
% load('chosenPoints20140728_Hela.mat')

%% Test transformation

im1t=imtransform(im_fix,tform,'XData', [1 size(im_mov,2)], 'YData', [1 size(im_mov,1)], ...
                'Size', size(im_mov));
figure;imshowpair(sqrt(double(im1t)),sqrt(double(im_mov)),'blend');axis equal;axis tight;

% show2subs(im1t,im_fix)

%% Choose region of interest
figure,
[~,rect]= imcrop(im1t); % rect = [x0, y0, widthx, widthy]
rect = round(rect);

%% Load stack of images to be transformed (1st stack) 
% file1 = 'C:\Users\Lukestom\Documents\ÈVUT\Ph.D. III.semestr\SOFI method\DATA\7\red\30msMEA20mMOSSGreenRedUV5mW_DU897_BV_1948.tif';
% file1 = 'K:\sharipov\20140529_Hela_3color\7\red\30msMEA20mMOSSGreenRedUV5mW_DU897_BV_1948.tif';
% Nframes = 100;
% imstack1 = load_tifFile(file1, Nframes);
% imstack1 = load_tifFile(file1);
% Nframes = size(imstack1,3);
% imstack1 = flipud(imstack1);%flip matrix up to down

% select tif to be loaded into imstack1 and imstack2

fn1 = 'data1.bin';
pn1 = 'K:\sharipov\20140729_Hela_from20140527\002\';
pn1 = 'H:\Sharipov\Measurements\20140728_Hela\002\';
fn2 = 'data2.bin';
pn2 = pn1;

Nframes =8000;
imstack1 = loadbinfile(pn1,fn1,Nframes);

% [fn1,pn1]=uigetfile('*.tif','MultiSelect','on');
% [fn2,pn2]=uigetfile('*.tif','MultiSelect','on');
% imstack1 = load_multitifFile(fn1,pn1);
%% Transform stack of images
Nframes = size(imstack1,3);
% imstack1t = zeros([size(imc2),Nframes]);

for ii = 1:Nframes

im_fix = imstack1(:,:,ii);    
im1t=imtransform(im_fix,tform,'XData', [1 size(im_mov,2)], 'YData', [1 size(im_mov,1)], ...
                'Size', size(im_mov));
            
imstack1t(:,:,ii) = im1t;
disp(['processing im:', num2str(ii)]);

end
% clear imstack1;

% crop image stack 
imstack1crop = imstack1t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);
clear imstack1t;

%% Load second stack of images
% file2 = 'C:\Users\Lukestom\Documents\ÈVUT\Ph.D. III.semestr\SOFI method\DATA\7\green\30msMEA20mMOSSGreenRedUV5mW_Luc247_MONO_0667.tif';
% file2 = 'K:\sharipov\20140529_Hela_3color\7\green\30msMEA20mMOSSGreenRedUV5mW_Luc247_MONO_0667.tif';

%%% loading tiff files and multitif files
% imstack2 = load_tifFile(file2, Nframes);
% imstack2 = flipud(imstack2);%flip matrix up to down
% imstack2 = load_multitifFile(fn2,pn2);
% imstack2 = flipdim(imstack2,1);%flip matrix up to down

imstack2 = loadbinfile(pn2,fn2,Nframes);
imstack2=flipdim(imstack2,2);
imstack2=flipdim(imstack2,1);

% crop image stack 
imstack2crop = imstack2(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);
clear imstack2;

% show2subs(std(double(imstack1crop),0,3),std(double(imstack2crop),0,3));
figure;imshowpair(sqrt(double(mean(imstack1crop,3))),sqrt(mean(imstack2crop,3)),'falsecolor');axis equal;axis tight;

%% Calculate MC SOFI cumulants

subseqlength=500;
totlength=Nframes;
start=501; % skip first n-images
Nss=floor((totlength-start+1)/subseqlength); % number of subsequences
orders=1:3;
for ns=1:Nss
    disp(['processing Nss:',num2str(ns)]);
    fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
    sofi=mcSofiCumulants(permute(cat(4,imstack1crop(:,:,fr),imstack2crop(:,:,fr)),[1 2 4 3]),[],[],[],orders);
    if ns==1
        c=sofi;
        c=cellfun(@(x)repmat(x,[1 1 1 Nss]),c,'UniformOutput',0);
    else
        for io=orders
            c{io}(:,:,:,ns)=sofi{io};
        end
    end
end

%% Calculate 2D cross-cumulants of non-transformed raw data (1st stack)
for ns=1:Nss
    disp(['processing Nss:',num2str(ns)]);
    fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
    sofi=mcSofiCumulantsSingle(imstack1(:,:,fr),[],[],[],orders);
    if ns==1
        c2=sofi;
        c2=cellfun(@(x)repmat(x,[1 1 Nss]),c2,'UniformOutput',0);
    else
        for io=orders
            c2{io}(:,:,ns)=sofi{io};
        end
    end
end

%% Flattening
cflat=c;
cm=cellfun(@(x)mean(x,4),c,'UniformOutput',0);
for ic=1:size(c{2},3)-1
    disp(['Flattening ic:',num2str(ic)]);
    for io=2:3
        for xx=0:io-1
            for yy=0:io-1
                cflat{io}(yy+1:io:end,xx+1:io:end,ic,:)=cflat{io}(yy+1:io:end,xx+1:io:end,ic,:)./repmat(std(reshape(cm{io}(yy+1:io:end,xx+1:io:end,ic),[],1)),size(cflat{io}(yy+1:io:end,xx+1:io:end,ic,:)));
            end
        end
        cflat{io}(:,:,ic,:)=cflat{io}(:,:,ic,:)*mean(mean(cm{io}(:,:,ic)))/mean(reshape(cflat{io}(:,:,ic,:),[],1));
    end
end
%% Flattening of 2D cross-cumulants of non-transformed raw data
c2flat=c2;
cm=cellfun(@(x)mean(x,3),c2,'UniformOutput',0);
for io=2:3
    for xx=0:io-1
        for yy=0:io-1
            c2flat{io}(yy+1:io:end,xx+1:io:end,:)=c2flat{io}(yy+1:io:end,xx+1:io:end,:)./repmat(std(reshape(cm{io}(yy+1:io:end,xx+1:io:end),[],1)),size(c2flat{io}(yy+1:io:end,xx+1:io:end,:)));
            
        end
    end
    c2flat{io}=c2flat{io}*mean(cm{io}(:))/mean(c2flat{io}(:));
end

%% Transform second channel cross-cumulants and combine cumulants
c2t=c2flat;
for io=orders
    T=tform.tdata.T;
    T(3,1:2)=io*T(3,1:2);
    T(3,2)=T(3,2)-io; %might have to be corrected differently - case by case
    c2t{io}=imtransform(c2flat{io},maketform('affine',T),'XData',[1 size(c{io},2)],'YData',[1 size(c{io},1)],'Size',[size(c{io},1) size(c{io},2)]);
end

tmp=cellfun(@(x)reshape(x(end:-1:1,:,:),size(x,1),size(x,2),1,[]),c2t,'UniformOutput',0);
for io=orders
    cflat{io}(:,:,end,:)=tmp{io};
end

%% Spectral unmixing - using 2nd order cumulant - 3 colors
% R1=0.25;%0.26;
% R2=0.9;%0.96;
% R3=1000;

R1=0.2;%0.26;
R2=0.01;%0.96;
R3=300;

R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=R^(-1);

mc2=cflat{2};
% nth=zeros(size(mc2,1),size(mc2,2),2);
umx1c2=squeeze(Rinv(1,1)*mc2(:,:,3,:)+Rinv(1,2)*mc2(:,:,2,:)+Rinv(1,3)*mc2(:,:,1,:));
umx2c2=squeeze(Rinv(2,1)*mc2(:,:,3,:)+Rinv(2,2)*mc2(:,:,2,:)+Rinv(2,3)*mc2(:,:,1,:));
umx3c2=squeeze(Rinv(3,1)*mc2(:,:,3,:)+Rinv(3,2)*mc2(:,:,2,:)+Rinv(3,3)*mc2(:,:,1,:));

% figure;imagesc(mean(umx1c2,3),[0 max(umx1c2(:))]);
% figure;imagesc(mean(umx2c2,3),[0 max(umx2c2(:))]);
% figure;imagesc(mean(umx3c2,3),[0 max(umx3c2(:))]);

% Merging unmixed cumulants into RGB image

col1=mean(umx1c2,3);
col2=mean(umx2c2,3);
col3=mean(umx3c2,3);
r=abs(col3).^(1/1);
g=abs(col2).^(1/1);
b=abs(col1).^(1/1);
r=r/max(r(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=b/max(b(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');


load myRmap;
load myGmap;
load myBmap;

% r = flipdim(r,2);
r = flipdim(r,1);

figure, 

subplot(131);imshow(b,[0 0.5]);colormap(myBmap);freezeColors;
subplot(132);imshow(g,[0 0.5]);colormap(myGmap);freezeColors;
subplot(133);imshow(r,[0 0.5]);colormap(myRmap);freezeColors;

mergedcum=cat(3,r,g,b);
mergedcum=1.5*mergedcum;
mergedcum(mergedcum>1)=1;
figure;
image(mergedcum);
% imwrite(uint8(2^8*merged),'merged_unmixed.png');
%%

ims1=mean(squeeze(mc2(:,:,2,:)),3);
ims1=ims1-min(ims1(:));
ims1=ims1/max(ims1(:));
% ims1 = flipud(ims1);
imwrite(ims1,'ims2.tif');

%% unmix 3 colors 3rd order

% R1=0.25;%0.26;
% R2=0.9;%0.96;
% R3=1000;

R1=0.2;%0.26;
R2=0.001;%0.96;
R3=300;

R=[R1^3 R2^3 R3^3; R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=pinv(R);

mc3=cflat{3};
nth=zeros(size(mc3,1),size(mc3,2),2);
umx1c3=squeeze(Rinv(1,1)*mc3(:,:,4,:)+Rinv(1,2)*mc3(:,:,3,:)+Rinv(1,3)*mc3(:,:,2,:)+Rinv(1,4)*mc3(:,:,1,:));
umx2c3=squeeze(Rinv(2,1)*mc3(:,:,4,:)+Rinv(2,2)*mc3(:,:,3,:)+Rinv(2,3)*mc3(:,:,2,:)+Rinv(2,4)*mc3(:,:,1,:));
umx3c3=squeeze(Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:)+Rinv(3,4)*mc3(:,:,1,:));

% R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
% Rinv=R^(-1);
% 
% mc3=c{3}; %mean(cflat{3},4)
% nth=zeros(size(mc3,1),size(mc3,2),2);
% umx1c3=squeeze(Rinv(1,1)*mc3(:,:,3,:)+Rinv(1,2)*mc3(:,:,2,:)+Rinv(1,3)*mc3(:,:,1,:));
% umx2c3=squeeze(Rinv(2,1)*mc3(:,:,3,:)+Rinv(2,2)*mc3(:,:,2,:)+Rinv(2,3)*mc3(:,:,1,:));
% 
% R=[R1^3 R2^3 R3^3; R1^2 R2^2 R3^2; R1 R2 R3];
% Rinv=R^(-1);
% umx3c3=squeeze(Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:));

col1=mean(umx1c3,3);
col2=mean(umx2c3,3);
col3=mean(umx3c3,3);
r=abs(col3).^(1/3);
g=abs(col2).^(1/3);
b=abs(col1).^(1/3);
r=r/max(r(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=b/max(b(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

merged=cat(3,r,g,b);
merged=2*merged;
merged(merged>1)=1;

figure, 

subplot(131);imshow(b,[0 0.5]);colormap(myBmap);freezeColors;
subplot(132);imshow(g,[0 0.5]);colormap(myGmap);freezeColors;
subplot(133);imshow(r,[0 0.5]);colormap(myRmap);freezeColors;

% figure;
% image(merged);

%% Linearization

[sy,sx,~] = size(umx1c3);
order=2;

% h2 = fspecial('Gaussian',[round(sy/2),round(sx/2)],2.3*sqrt(order)/2/sqrt(2*log(2))*sy);
h2 = fspecial('Gaussian',[sy,sx],2.3*sqrt(order)/2/sqrt(2*log(2)));
% h2 = padarray(h2,[sy-round(sy/2),sx-round(sx/2)]);
h = fspecial('Gaussian',29,2.3*sqrt(order)/2/sqrt(2*log(2)));



im1=0;
im2=0;
im3=0;

settings.gamma = 10000;
settings.beta = settings.gamma;
settings.alpha = 1;      
settings.reltol = 1e-4;
settings.maxIter = 50;
settings.Lp = 1;

for n=1:size(umx1c2,3)
    disp(['Linearization un:',num2str(n)]);
    
    tmp=abs(umx1c3(:,:,n));
%     tmp=medfilt2(tmp);
    scale=max(tmp(:));
%     tmp=(scale*deconvlucy(tmp/scale,fspecial('Gaussian',29,2.3*sqrt(order)/2/sqrt(2*log(2))),5)).^(1/order);
    
    [imOut,relChangeAll] = deconvAugLag3D(tmp/scale,fftn(h2,[sy,sx]),settings);
    tmp = scale*imOut;
    im1=im1+tmp;
    
    tmp=abs(umx2c3(:,:,n));
%     tmp=medfilt2(tmp);
    scale=max(tmp(:));
%     tmp=(scale*deconvlucy(tmp/scale,fspecial('Gaussian',29,2.5*sqrt(order)/2/sqrt(2*log(2))),5)).^(1/order);
    [imOut,relChangeAll] = deconvAugLag3D(tmp/scale,fftn(h2,[sy,sx]),settings);
    tmp = scale*imOut;
    im2=im2+tmp;
    
    tmp=abs(umx3c3(:,:,n));
%     tmp=medfilt2(tmp);
    scale=max(tmp(:));
%     tmp=(scale*deconvlucy(tmp/scale,fspecial('Gaussian',29,2.7*sqrt(order)/2/sqrt(2*log(2))),5)).^(1/order);
    [imOut,relChangeAll] = deconvAugLag3D(tmp/scale,fftn(h2,[sy,sx]),settings);
    tmp = scale*imOut;
    im3=im3+tmp;
end

% save('linumxcum.mat','im1','im2','im3');
% save('linumxcum2_transfCorrected_auglag.mat','im1','im2','im3');
%% Merging linearized deconvolved unmixed cumulants into final RGB image

r=abs(im3);
r=r-2*median(r(:));
r(r<0)=0;
g=abs(im2);
g=g-2*median(g(:));
g(g<0)=0;
b=abs(im1);
b=b-2*median(b(:));
b(b<0)=0;
r=r/max(r(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=b/max(b(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

merged=cat(3,r,g,b);
merged=1.5*merged;
merged(merged>1)=1;
% merged2 = ifftshift(merged);

figure;
image(merged);

show2subs(mergedcum,merged)

%%
r = im3;
% r = flipdim(r,2);
r = flipdim(r,1);
g = im2;
b = im1;
r=6*r/max(r(:));
g=9*g/max(g(:));
b=6*b/max(b(:));

r = r-1*median(r(:));
r(r<0)=0;
g = g-1*median(g(:));
g(g<0)=0;
b = b-1*median(b(:));
b(b<0)=0;

% r = medfilt2(r,[5 5]);
% g = medfilt2(g,[5 5]);
% b = medfilt2(b,[5 5]);

load myRmap;
load myGmap;
load myBmap;

figure, 

subplot(131);imshow(b,[0 0.5]);colormap(myBmap);freezeColors;
subplot(132);imshow(g,[0 0.5]);colormap(myGmap);freezeColors;
subplot(133);imshow(r,[0 0.5]);colormap(myRmap);freezeColors;

imrgb = cat(3,r,g,b);
figure, 
imshow(imrgb,[0 0.5]);

%% Save unmixed cumulants

% save('unmixedCumulants_2ord_transfCorrected_Hela3Colors','umx1c2','umx2c2','umx3c2')

%% Export to tif files

% nameout = 'imstack1t';
% 
% % widefield stack
% imwrite(uint16(round(imstack(:,:,1)*65536)), [nameout,'.tif'],'Compression','none')
% 
% for k = 2:size(imstack1t,3)
%     imwrite(uint16(round(imstack1t(:,:,k)*65536)), [nameout,'.tif'], 'Compression','none','writemode', 'append');
% end

%% Export to binary files

% fid = fopen('data1.bin', 'w');
% fwrite(fid, imstack, 'uint16');
% fclose(fid);





