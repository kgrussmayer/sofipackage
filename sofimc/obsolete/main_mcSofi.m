%%% run multicolor sofi
%%% adrien.descloux@epfl.ch 11/2018
%%% tomas.lukes@epfl.ch
%%% kristin.grussmayer@epfl.ch

% close all; clear; 

%
%%% Initialization
%

addpath(genpath('..\utils'));
addpath(genpath('funcs'));
config_MCSOFI2D;

%%% Load the calibration files
imc1 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc1, settings.io.fileExtension]);
imc2 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc2, settings.io.fileExtension,]);
if strcmp(settings.io.dataType, 'exp')
    imc2 = flip(imc2,2);
end

% imc1(imc1==0)=min(imc1(imc1>0));
% imc2(imc2==0)=min(imc2(imc2>0));

% what is indcom used for?
[~,indcom]=max(squeeze(mean(max(imc1))).*squeeze(mean(max(imc2))));
im_fix=imc1(:,:,indcom);
im_mov=imc2(:,:,indcom);

show2subs(im_fix,im_mov)

%
%%% Run calibration (register images from 2 color channels)
%

cal = runcalib_mcsofi(imc1,imc2,settings.sys,settings.cal,settings.io); 

%% Load stack of images to be transformed (1st stack) 
% settings.io.fileExtension = '.tif';
imstack1 = load_tifFile([settings.io.pn, filesep, settings.io.fn1, settings.io.fileExtension],settings.sys.sub);

imstack2 = load_tifFile([settings.io.pn, filesep, settings.io.fn2, settings.io.fileExtension],settings.sys.sub);
if strcmp(settings.io.dataType, 'exp')
    imstack2=flip(imstack2,2);
end

st1 = std(single(imstack1),[],3);
st2 = std(single(imstack2),[],3);

% select tiff to be loaded into imstack1 and imstack2
% [fn1,pn1]=uigetfile('*.tif','MultiSelect','on');
% [fn2,pn2]=uigetfile('*.tif','MultiSelect','on');
% imstack1 = load_multitifFile(fn1,pn1);

%% Region of interest

if isfield(cal,'rect')  
    rect = cal.rect;
else
    rect = [1,1,size(imstack1,2),size(imstack1,1)]; %like imcrop, x,y,width,height
end

%% Transform stack of images
Nframes = size(imstack1,3);
% imstack1t = zeros([size(imc2),Nframes]);
imstack2t = zeros(size(imstack1),'single');

% if the data are exp, need to compute the proper coregistration from data
if strcmp(settings.io.dataType, 'exp')
    st2t = imtranslate(st2,[0 -cal.tf{1}.B(1)]);
    temp = zeros(size(st1,1),2*size(st1,2)-1);
    for k = 1:size(st1,1)
        temp(k,:) = normxcorr(st1(k,:),st2t(k,:));
    end
    temp(isnan(temp)) = 0;
    temp = mean(temp,1);
    [p,ind] = findpeaks(double(temp),'MinPeakWidth',2);
    [~,v] = max(abs(ind-size(st1,2)).*(abs(ind-size(st1,2))<settings.cal.maxshiftx));
    shiftx = ind(v)-size(st1,2);
    st2t = imtranslate(st2t,[-shiftx 0]);
end

% apply transform to data
for ii = 1:Nframes

im_mov = imstack2(:,:,ii);    
% im1t=imtransform(im_fix,tform,'XData', [1 size(im_mov,2)], 'YData', [1 size(im_mov,1)], ...
%                 'Size', size(im_mov));
% im1t = imwarp(im_mov, cal.tf{1},'OutputView',imref2d(size(im_mov)));
if strcmp(settings.io.dataType,'exp')
    im1t = imtranslate(im_mov,[-shiftx -cal.tf{1}.B(1)]);
else
    im1t = imtranslate(im_mov,[-cal.tf{1}.A(1) -cal.tf{1}.B(1)]); % hard coded values because calibration file is wrong!
end
imstack2t(:,:,ii) = single(im1t);
% disp(['processing im:', num2str(ii)]);

end

% compute coregistration mask and proper crop
mask = not(im1t(:,:,1) == 0);
[r,c] = find(mask>0);
rect = [min(c) min(r) max(c)-min(c) max(r)-min(r)];

% crop image stack 
imstack1crop = imstack1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);

% crop image stack 
imstack2crop = single(imstack2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:));

% show2subs(std(double(imstack1crop),0,3),std(double(imstack2crop),0,3));
figure;imshowpair(st1,st2t,'falsecolor');axis equal;axis tight;

%% Calculate MC SOFI cumulants
orders = settings.sys.orders;
subseqlength = settings.sys.wsize;
start = settings.sys.start;

Nss=floor((Nframes-start+1)/subseqlength); % number of subsequences

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
    sofi=mcSofiCumulantsSingle(imstack2(:,:,fr),[],[],[],orders);
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
    disp(io)
    
%     T=tform.tdata.T;
%     T(3,1:2)=io*T(3,1:2);
%     T(3,2)=T(3,2)-io; %might have to be corrected differently - case by case
%     c2t{io}=imtransform(c2flat{io},maketform('affine',T),'XData',[1 size(c{io},2)],'YData',[1 size(c{io},1)],'Size',[size(c{io},1) size(c{io},2)]);
%     c2t{io} = imwarp(c2flat{io}, cal.tf{io},'OutputView',imref2d(size(c2t{io})));
    c2t{io} = imtranslate(c2flat{io},io.*[-shiftx -cal.tf{1}.B(1)]);
    % crop image stack 
    if strcmp(settings.io.dataType,'exp')
     % For Real DATA
     c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-2*io:io*rect(2) + io*rect(4)-ceil(3*io/2)-2*io,...
         io*rect(1)+floor(3*io/2)-2*io:io*rect(1)+io*rect(3)-ceil(3*io/2)-2*io,:);
    else
    % For SIMULATION
    c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-2*io:io*rect(2) + io*rect(4)-ceil(3*io/2)-io,...
        io*rect(1)+floor(3*io/2)-io:io*rect(1)+io*rect(3)-ceil(3*io/2)-io,:);
    end
end

tmp=cellfun(@(x)reshape(x(:,:,:),size(x,1),size(x,2),1,[]),c2t,'UniformOutput',0);
for io=orders
    cflat{io}(:,:,end,:)=tmp{io};
end

%% Spectral unmixing - using 2nd order cumulant - 3 colors
R1 = settings.mc.R1;
R2 = settings.mc.R2;
R3 = settings.mc.R3;

R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=R^(-1);

mc2=cflat{2};
% mc2=cflat{2};
% nth=zeros(size(mc2,1),size(mc2,2),2);
umx1c2=squeeze(Rinv(1,1)*mc2(:,:,3,:)+Rinv(1,2)*mc2(:,:,2,:)+Rinv(1,3)*mc2(:,:,1,:));
umx2c2=squeeze(Rinv(2,1)*mc2(:,:,3,:)+Rinv(2,2)*mc2(:,:,2,:)+Rinv(2,3)*mc2(:,:,1,:));
umx3c2=squeeze(Rinv(3,1)*mc2(:,:,3,:)+Rinv(3,2)*mc2(:,:,2,:)+Rinv(3,3)*mc2(:,:,1,:));
% umx3c2 = flip(umx3c2,1);

% figure;imagesc(mean(umx1c2,3),[0 max(umx1c2(:))]);
% figure;imagesc(mean(umx2c2,3),[0 max(umx2c2(:))]);
% figure;imagesc(mean(umx3c2,3),[0 max(umx3c2(:))]);

% figure;imagesc(mean(mc2(:,:,1,1),3));
% figure;imagesc(mean(mc2(:,:,2,1),3));
% figure;imagesc(mean(mc2(:,:,3,1),3));

% Merging unmixed cumulants into RGB image
lincoeff = 1;%/0.6;
col1=mean(umx1c2,3);
col2=mean(umx2c2,3);
col3=mean(umx3c2,3);
r=(col3(:,:));%.^lincoeff;
g=(col2);%.^lincoeff;
b=(col1);%.^lincoeff;
r=r/max(r(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=b/max(b(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

% load myRmap;
% load myGmap;
% load myBmap;
myRmap = zeros(64,3); myRmap(:,1) = linspace(0,1,64); myRmap(:,3) = linspace(0,1,64);
myGmap = zeros(64,3); myGmap(:,2) = linspace(0,1,64); myGmap(:,1) = linspace(0,1,64);
myBmap = zeros(64,3); myBmap(:,3) = linspace(0,1,64);myBmap(:,2) = linspace(0,1,64);

% r = flipdim(r,2);
% r = flipdim(r,1);

% figure('units','normalized','outerposition',[0 0 1 1]);
figure(1235)
subplot(131);imshow(b,[0 1]);colormap(myBmap);freezeColors;title('channel 1')
subplot(132);imshow(g,[0 1]);colormap(myGmap);freezeColors;title('virtual channel') 
subplot(133);imshow(r,[0 1]);colormap(myRmap);freezeColors;title('channel 2')

% if settings.io.figsave == 1 
%             saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
%                 [settings.io.fn1,'_results_rgb_sofiraw2'],settings.io.figformat);
% end
% 
mergedcum=cat(3,r,g,b);
mergedcum=1.5*mergedcum;
mergedcum(mergedcum>1)=1;
figure;
image(mergedcum);
% imwrite(uint8(2^8*merged),'merged_unmixed.png');
%%

% ims1=mean(squeeze(mc2(:,:,2,:)),3);
% ims1=ims1-min(ims1(:));
% ims1=ims1/max(ims1(:));
% % ims1 = flipud(ims1);
% imwrite(ims1,'ims2.tif');

%% unmix 3 colors 3rd order

R=[R1^3 R2^3 R3^3; R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=pinv(R);

mc3=cflat{3};
nth=zeros(size(mc3,1),size(mc3,2),2);
umx1c3=squeeze(Rinv(1,1)*mc3(:,:,4,:)+Rinv(1,2)*mc3(:,:,3,:)+Rinv(1,3)*mc3(:,:,2,:)+Rinv(1,4)*mc3(:,:,1,:));
umx2c3=squeeze(Rinv(2,1)*mc3(:,:,4,:)+Rinv(2,2)*mc3(:,:,3,:)+Rinv(2,3)*mc3(:,:,2,:)+Rinv(2,4)*mc3(:,:,1,:));
umx3c3=squeeze(Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:)+Rinv(3,4)*mc3(:,:,1,:));
% umx3c3 = flip(umx3c3,1);
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
r=(col3).^(1/3);
g=(col2).^(1/3);
b=(col1).^(1/3);
r=r/max(r(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=b/max(b(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

merged=cat(3,r,g,b);
merged=1.5*merged;
merged(merged>1)=1;

figure, 

subplot(131);imshow(b,[0 1]);colormap(myBmap);freezeColors;
subplot(132);imshow(g,[0 1]);colormap(myGmap);freezeColors;
subplot(133);imshow(r,[0 1]);colormap(myRmap);freezeColors;

% figure;
% image(merged);

%% unmix 4 colors 3rd order
% missing, lost? 
% from obsolete folder
R=[R1^3 R2^3 R3^3 R4^3; R1^2 R2^2 R3^2 R4^2; R1 R2 R3 R4; 1 1 1 1];
Rinv=pinv(R);

% % for 3 color unmixing
% mc3=cflat{3};
% nth=zeros(size(mc3,1),size(mc3,2),2);
% umx1c3=squeeze(Rinv(1,1)*mc3(:,:,4,:)+Rinv(1,2)*mc3(:,:,3,:)+Rinv(1,3)*mc3(:,:,2,:)+Rinv(1,4)*mc3(:,:,1,:));
% umx2c3=squeeze(Rinv(2,1)*mc3(:,:,4,:)+Rinv(2,2)*mc3(:,:,3,:)+Rinv(2,3)*mc3(:,:,2,:)+Rinv(2,4)*mc3(:,:,1,:));
% umx3c3=squeeze(Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:)+Rinv(3,4)*mc3(:,:,1,:));

mc3=cflat{3};
nth=zeros(size(mc3,1),size(mc3,2),2);
%unmixing 4 colors using 3rd order, unmixing from obsolete folder,
%variables exchanged according to 3 color nomenclature above
umx1c3=Rinv(1,1)*mc3(:,:,4,:)+Rinv(1,2)*mc3(:,:,3,:)+Rinv(1,3)*mc3(:,:,2,:)+Rinv(1,4)*mc3(:,:,1,:);
umx2c3=Rinv(2,1)*mc3(:,:,4,:)+Rinv(2,2)*mc3(:,:,3,:)+Rinv(2,3)*mc3(:,:,2,:)+Rinv(2,4)*mc3(:,:,1,:);
umx3c3=Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:)+Rinv(3,4)*mc3(:,:,1,:);
umx4c3=Rinv(4,1)*mc3(:,:,4,:)+Rinv(4,2)*mc3(:,:,3,:)+Rinv(4,3)*mc3(:,:,2,:)+Rinv(4,4)*mc3(:,:,1,:);

% further processing is still missing

%% Linearization

[sy,sx,~] = size(umx1c2);
order=2;

% h2 = fspecial('Gaussian',[round(sy/2),round(sx/2)],2.3*sqrt(order)/2/sqrt(2*log(2))*sy);
h2 = fspecial('Gaussian',[sy,sx],2.3*sqrt(order)/2/sqrt(2*log(2)));
% h2 = padarray(h2,[sy-round(sy/2),sx-round(sx/2)]);
h = fspecial('Gaussian',29,2.3*sqrt(order)/2/sqrt(2*log(2)));

settings.gamma = 2000;
settings.beta = settings.gamma;
settings.alpha = 1;      
settings.reltol = 1e-4;
settings.maxIter = 20;
settings.Lp = 1;

for n=1:size(umx1c2,3)
    disp(['Linearization un:',num2str(n)]);
    
    tmp=(umx1c2(:,:,n));
    tmp(tmp<0)=0;
%     tmp=medfilt2(tmp);
    scale=max(tmp(:));
%     tmp=(scale*deconvlucy(tmp/scale,fspecial('Gaussian',29,2.3*sqrt(order)/2/sqrt(2*log(2))),5)).^(1/order);
    
    [imOut,relChangeAll] = deconvAugLag3D(tmp/scale,fftn(h2,[sy,sx]),settings);
    tmp = scale*imOut;
    im1(:,:,n)=tmp;
    
    tmp=(umx2c2(:,:,n));
    tmp(tmp<0)=0;
%     tmp=medfilt2(tmp);
    scale=max(tmp(:));
%     tmp=(scale*deconvlucy(tmp/scale,fspecial('Gaussian',29,2.5*sqrt(order)/2/sqrt(2*log(2))),5)).^(1/order);
    [imOut,relChangeAll] = deconvAugLag3D(tmp/scale,fftn(h2,[sy,sx]),settings);
    tmp = scale*imOut;
    im2(:,:,n)=tmp;
    
    tmp=(umx3c2(:,:,n));
    tmp(tmp<0)=0;
%     tmp=medfilt2(tmp);
    scale=max(tmp(:));
%     tmp=(scale*deconvlucy(tmp/scale,fspecial('Gaussian',29,2.7*sqrt(order)/2/sqrt(2*log(2))),5)).^(1/order);
    [imOut,relChangeAll] = deconvAugLag3D(tmp/scale,fftn(h2,[sy,sx]),settings);
    tmp = scale*imOut;
    im3(:,:,n)=tmp;
end
%%
figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(131);imshow(mean(im1,3),[]);colormap(myBmap);freezeColors;title('channel 1')
subplot(132);imshow(mean(im2,3),[]);colormap(myGmap);freezeColors;title('virtual channel') 
subplot(133);imshow(mean(im3,3),[]);colormap(myRmap);freezeColors;title('channel 2')

if settings.io.figsave ==1 
            saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
                [settings.io.fn1,'_results_rgb_sofidec2'],settings.io.figformat);
end
% save('linumxcum.mat','im1','im2','im3');
% save('linumxcum2_transfCorrected_auglag.mat','im1','im2','im3');
%% Merging linearized deconvolved unmixed cumulants into final RGB image
bcgsub = 0.5;
r=mean(im3,3).^lincoeff;
% r=r-bcgsub*median(r(:));
% r(r<0)=0;
g=mean(im2,3).^lincoeff;
% g=g-bcgsub*median(g(:));
% g(g<0)=0;
b=mean(im1,3).^lincoeff;
% b=b-bcgsub*median(b(:));
% b(b<0)=0;
r=r/max(r(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=b/max(b(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');


merged=cat(3,r,g,b);
merged=1.5*merged;
% merged(merged>1)=1;
% merged2 = ifftshift(merged);
figure(902345);image(merged)
figure(340530);
subplot(131);imagesc(r); title('Red')
subplot(132);imagesc(g); title('Green')
subplot(133);imagesc(b); title('Blue')

%%
% figure;
% image(merged);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(131);
imshowpair(((mean(imstack1crop,3))),(mean(imstack2crop,3)),'falsecolor');axis equal;axis tight;
title('Averaged images, the two input channels overlaid')
subplot(132);
imshow(mergedcum,[])
title('Spectral-cross cumulants 2nd order')
subplot(133);
imshow(merged,[]);
title('Spectral-cross cumulants 2nd order, deconvolved, linearized')

if settings.io.figsave ==1 
            saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
                [settings.io.fn1,'_results_rgb_composite_avg_sofiraw2_sofilin2'],settings.io.figformat);
end

% show2subs(mergedcum,merged)

%%
% r = im3;
% % r = flipdim(r,2);
% r = flipdim(r,1);
% g = im2;
% b = im1;
% r=6*r/max(r(:));
% g=9*g/max(g(:));
% b=6*b/max(b(:));
% 
% r = r-1*median(r(:));
% r(r<0)=0;
% g = g-1*median(g(:));
% g(g<0)=0;
% b = b-1*median(b(:));
% b(b<0)=0;
% 
% % r = medfilt2(r,[5 5]);
% % g = medfilt2(g,[5 5]);
% % b = medfilt2(b,[5 5]);
% 
% load myRmap;
% load myGmap;
% load myBmap;
% 
% figure, 
% 
% subplot(131);imshow(b,[0 0.5]);colormap(myBmap);freezeColors;
% subplot(132);imshow(g,[0 0.5]);colormap(myGmap);freezeColors;
% subplot(133);imshow(r,[0 0.5]);colormap(myRmap);freezeColors;
% 
% imrgb = cat(3,r,g,b);
% figure, 
% imshow(imrgb,[0 0.5]);

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





