%%% test code for optimizing MCSOFI processing on specific data
%%% adrien.descloux@epfl.ch 18/06/2018

close all; clear; 

%
%%% Initialization

addpath(genpath('..\utils'));
addpath(genpath('funcs'));
config_mc;

% select folder to batch process

pname = uigetdir('K:\Kristin');
% pname = 'K:\Kristin\20180201'; % to be removed 
settings.io.pn = pname;

list = dir(pname);
list(1:2) = []; % remove . and .. path from the list

fID = '013'; % specify which file is under investigation

% look for calibration file
isCal = 0;
for k = length(list):-1:1 % loop backwards because we are changing the list on the fly
    if strfind(list(k).name,'calibration')
        settings.io.pnc = [pname,filesep,list(k).name];
        isCal = 1;
        list(k) = []; % remove calibration folder from the list
    elseif strfind(list(k).name,'results')
        list(k) = []; % remove results folder from the list
    elseif isempty(strfind(list(k).name,'.tiff'))&& isempty(strfind(list(k).name,'.tif'))
        list(k) = []; % remove extra folders from tomas procesing
    elseif strfind(list(k).name,'.mat')
        list(k) = [];
    end
end
if isCal == 0
    disp('Calibration file not found !');
    disp('Please make sure the folder is properly organized!'); 
else

% make list of paired files
ex = length(list(1).name)-strfind(list(1).name,'.');
idMax = 0; listC = [];
for k = 1:length(list)
    temp = str2double(list(k).name(end-ex-3:end-ex-1));
    if temp > idMax
        idMax = temp;
    end
end

rList = []; tList = [];
for k = 0:idMax 
    id = num2str(k); while length(id) < 3; id = ['0',id];end
    id = ['_',id];
    for h = length(list):-1:1
       if strfind(list(h).name,id) % if we find proper ID number
           if strfind(list(h).name,'reflected')
               rList(end+1).name = list(h).name;
               rList(end).id = id(2:end);
           elseif strfind(list(h).name,'transmitted')
               tList(end+1).name = list(h).name;
               tList(end).id = id(2:end);
           end
       end
    end
end

if length(tList) ~= length(rList)
    warning('Error number of reflected and transmitted files are different');
end
    rFold = [];
    fpos = [];
   for k = 1:length(rList)
       if strcmp(rList(k).id,fID)
           % make new results folder
            rFold = ['results_test_MCSOFI_',fID,'_',getID];
            mkdir(pname,rFold);
            fpos = k;
            break;
       end
   end
    if isempty(rFold)
        disp(['Unable to find specified file ID :',fID])
    else
    % summary of folder to be processed
        disp(['Selected folder : ',settings.io.pn])
        disp(['Calibration folder : ',settings.io.pnc])
        disp(['Number of files detected : ',num2str(length(rList))])
    end
end

%% DATA loading and coregistration 

k = fpos;
    disp(['Processing file # ',num2str(k),', ID : ',rList(k).id])
    
%%% Load the calibration files
imc1 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc1, settings.io.fileExtension]);
imc2 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc2, settings.io.fileExtension]);
if strcmp(settings.io.dataType, 'exp')
    imc2 = flip(imc2,2);
end

% imc1(imc1==0)=min(imc1(imc1>0));
% imc2(imc2==0)=min(imc2(imc2>0));

% what is indcom used for?
[~,indcom]=max(squeeze(mean(max(imc1))).*squeeze(mean(max(imc2))));
im_fix=imc1(:,:,indcom);
im_mov=imc2(:,:,indcom);

show2subs(im_fix,im_mov,1)
close(1)

%%% Run calibration (register images from 2 color channels)
disp('Compute calibration')
settings.cal.roix = []; settings.cal.roiy = [];
cal = runcalib_mcsofi(imc1,imc2,settings.sys,settings.cal,settings.io); 

%% Load stack of images to be transformed (1st stack) 
% settings.io.fileExtension = '.tif';
settings.io.fn1 = rList(k).name;
settings.io.fn2 = tList(k).name;
disp('Loading reflected data')
imstack1 = load_tifFile([settings.io.pn, filesep, settings.io.fn1],settings.sys.sub);
disp('Loading transmitted data')
imstack2 = load_tifFile([settings.io.pn, filesep, settings.io.fn2],settings.sys.sub);
if strcmp(settings.io.dataType, 'exp')
    imstack2=flip(imstack2,2);
end

disp('Compute stacks STD')
% compute the std of the stacks
st1 = std(single(imstack1(:,:,round(0.2*size(imstack1,3)):round(0.8*size(imstack1,3)))),[],3);
st2 = std(single(imstack2(:,:,round(0.2*size(imstack1,3)):round(0.8*size(imstack1,3)))),[],3);

%% Transform stack of images
disp('Compute correct coregistration')
Nframes = size(imstack1,3);
% imstack1t = zeros([size(imc2),Nframes]);
imstack2t = zeros(size(imstack1),'single');
% 2D high pass filter
kernel1 = ones(3)/9;
t1 = st1(50:end-50,:); t2 = st2(50:end-50,:); 
t1 = medfilt2(t1,[2 2]); t2 = medfilt2(t2,[2 2]); % remove salt&pepper noise
t1 = t1-imgaussfilt(single(t1),5); t2 = t2-imgaussfilt(single(t2),5); % remove background before cross-correlation

% st2t = imtranslate(st2, [0 -cal.tf{1}.B(1)]);
% t2 = imtranslate(t2, [0 -cal.tf{1}.B(1)]);
if strcmp(settings.io.dataType, 'exp')

    temp = xcorr2(t1,t2);
    temp((size(temp,1)+1)/2,:) = 0;
    temp(:,(size(temp,2)+1)/2) = 0;
    [~,ind] = max(temp(:));
    [y,x] = ind2sub(size(temp),ind);
%     temp(size(st1,1)) = 0; temp(1) = 0; temp(end) = 0;
%     [~,ind] = max(double(temp));
%     [~,v] = max(abs(ind-size(st1,2)).*(abs(ind-size(st1,2))<settings.cal.maxshiftx));
    shiftx = -(x-(size(temp,2)+1)/2);
    shifty = -(y-(size(temp,1)+1)/2);
    st2t = imtranslate(st2,[-shiftx -shifty]);
else
    shiftx = 0;
    shifty = 0;
    st2t = st2;
end

    
disp('Transform of imstack2')
% apply transform to data
for ii = 1:Nframes

im_mov = imstack2(:,:,ii);    
% im1t=imtransform(im_fix,tform,'XData', [1 size(im_mov,2)], 'YData', [1 size(im_mov,1)], ...
%                 'Size', size(im_mov));
% im1t = imwarp(im_mov, cal.tf{1},'OutputView',imref2d(size(im_mov)));
% if strcmp(settings.io.dataType,'exp')
%     im1t = imwarp(im_mov,cal.tf{1},'OutputView',imref2d(size(im_mov)));
    im1t = imtranslate(im_mov,[-shiftx -shifty]);
% else
%     im1t = imtranslate(im_mov,[-cal.tf{1}.A(1) -cal.tf{1}.B(1)]); % use normal coregistration value for 'sim'
% end
imstack2t(:,:,ii) = single(im1t);
end
settings.cal.shiftx = -shiftx;
settings.cal.shifty = -shifty; % -cal.tf{1}.B(1)

% compute coregistration mask and proper crop
mask = not(im1t(:,:,1) == 0);
[r,c] = find(mask>0);
rect = [min(c)+1 min(r)+1 max(c)-min(c)-2 max(r)-min(r)-2];
% SOFI processing requires even dimensions or you gets weird cropping
if mod(rect(3),2) == 0; rect(3) = rect(3)-1; end
if mod(rect(4),2) == 0; rect(4) = rect(4)-1; end

settings.cal.roix = [rect(1) rect(3)];
settings.cal.roiy = [rect(2) rect(4)];

% crop image stack 
imstack1crop = imstack1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);
% crop image stack 
imstack2crop = single(imstack2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:));

st1c = st1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3));
st2tc = st2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3)); fake = [];
fake(:,:,2) = 255.*st1c./max(st1c(:)); fake(:,:,3) = 255.*st2tc./max(st2tc(:)); fake(:,:,1) = 255.*st2tc./max(st2tc(:));
figure(2);imagesc(uint8(fake))
% figure(2);imshowpair(st1,st2t,'falsecolor');axis equal;axis tight;

%%

start = 1201;
stop = size(imstack1,3) - 500;
wsize = 5000; 150:200:950;

for j = 1:length(wsize)

% compute a unique tag to avoid overwriting saved data
tag = num2str(floor(length(dir([pname,filesep,rFold]))/14)); while length(tag) < 3; tag = ['0',tag];end
disp(['Processing tag :',tag])
% modify the any settings parameters here
settings.sys.start = start;
settings.sys.stop = stop;
settings.sys.wsize = wsize(j);

% The settings structure is automatically saved in a txt file with the
% processed results !!!


%% Calculate MC SOFI cumulants
disp('Compute SOFI cumulants')
orders = settings.sys.orders;
subseqlength = settings.sys.wsize;
start = settings.sys.start;
stop = settings.sys.stop;

Nss=floor((stop-start+1)/subseqlength); % number of subsequences
settings.sys.Nss = Nss;

for ns=1:Nss
%     disp(['processing Nss:',num2str(ns)]);
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
disp('Compute SOFI cumulant of imstack2')
for ns=1:Nss
%     disp(['processing Nss:',num2str(ns)]);
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
disp('Cumulant flattening')
cflat=c;
cm=cellfun(@(x)mean(x,4),c,'UniformOutput',0);
for ic=1:size(c{2},3)-1
%     disp(['Flattening ic:',num2str(ic)]);
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
disp('Transform and crop second SOFI channel')
c2t=c2flat;
off = [1 4 6];
for io=orders
%     disp(io)
%     T=tform.tdata.T;
%     T(3,1:2)=io*T(3,1:2);
%     T(3,2)=T(3,2)-io; %might have to be corrected differently - case by case
%     c2t{io}=imtransform(c2flat{io},maketform('affine',T),'XData',[1 size(c{io},2)],'YData',[1 size(c{io},1)],'Size',[size(c{io},1) size(c{io},2)]);
%      c2t{io} = imwarp(c2flat{io}, cal.tf{io},'OutputView',imref2d(size(c2t{io})));
    c2t{io} = imtranslate(c2flat{io},io.*[-shiftx -shifty]); % -cal.tf{1}.B(1)
    % crop image stack 
%     if strcmp(settings.io.dataType,'exp')
     % For Real DATA
     c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-off(io):io*rect(2) + io*rect(4)-ceil(3*io/2)-off(io),...
         io*rect(1)+floor(3*io/2)-off(io):io*rect(1)+io*rect(3)-ceil(3*io/2)-off(io),:);
%     else
%     % For SIMULATION
%         c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-2*io:io*rect(2) + io*rect(4)-ceil(3*io/2)-io,...
%             io*rect(1)+floor(3*io/2)-io:io*rect(1)+io*rect(3)-ceil(3*io/2)-io,:);
%     end
end

disp('Combine cumulants')
tmp=cellfun(@(x)reshape(x(:,:,:),size(x,1),size(x,2),1,[]),c2t,'UniformOutput',0);
for io=orders
    cflat{io}(:,:,end,:)=tmp{io};
end

%% Spectral unmixing - using 2nd order cumulant - 3 colors
disp('Spectral unmixing')
% R1 = settings.mc.R1;
% R2 = settings.mc.R2;
% R3 = settings.mc.R3;

% get R1,R2,R3 from calibration and filename information
if strcmp(settings.io.dataType,'exp')
[R1,R2,R3] = getUnmixCoef(settings.io.pnc,settings.io.fn1);
else
    load([settings.io.pn,filesep,'unmixCoef_',rList(k).id]) ; % load R1 R2 and R3 into the workspace
end
settings.mc.R1 = R1; settings.mc.R2 = R2; settings.mc.R3 = R3;

R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=R^(-1);

mc2=cflat{2};
% nth=zeros(size(mc2,1),size(mc2,2),2);
umx1c2=squeeze(Rinv(1,1)*mc2(:,:,3,:)+Rinv(1,2)*mc2(:,:,2,:)+Rinv(1,3)*mc2(:,:,1,:));
umx2c2=squeeze(Rinv(2,1)*mc2(:,:,3,:)+Rinv(2,2)*mc2(:,:,2,:)+Rinv(2,3)*mc2(:,:,1,:));
umx3c2=squeeze(Rinv(3,1)*mc2(:,:,3,:)+Rinv(3,2)*mc2(:,:,2,:)+Rinv(3,3)*mc2(:,:,1,:));


% Merging unmixed cumulants into RGB image
lincoeff = 1;
settings.mc.lincoeff2 = lincoeff;
col1=mean(umx1c2,3); col1(col1<0) = 0;
col2=mean(umx2c2,3); col2(col2<0) = 0;
col3=mean(umx3c2,3); col3(col3<0) = 0;

r=(col3).^lincoeff;
g=(col2).^lincoeff;
b=(col1).^lincoeff;
r=r/max(r(:));
g=g/max(g(:));
b=b/max(b(:));


myRmap = zeros(64,3); myRmap(:,1) = linspace(0,1,64); myRmap(:,3) = linspace(0,1,64);
myGmap = zeros(64,3); myGmap(:,2) = linspace(0,1,64); myGmap(:,1) = linspace(0,1,64);
myBmap = zeros(64,3); myBmap(:,3) = linspace(0,1,64);myBmap(:,2) = linspace(0,1,64);


% figure(4)
% subplot(131);imshow(b,[0 1]);colormap(myBmap);freezeColors;title('channel 1')
% subplot(132);imshow(g,[0 1]);colormap(myGmap);freezeColors;title('virtual channel') 
% subplot(133);imshow(r,[0 1]);colormap(myRmap);freezeColors;title('channel 2')


mergedcum=cat(3,r,g,b);
mergedcum=1.5*mergedcum;
mergedcum(mergedcum>1)=1;
figure(5);
image(mergedcum);
set(gca,'xtick',[],'ytick',[])

%% unmix 3 colors 3rd order

R=[R1^3 R2^3 R3^3; R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=pinv(R);

mc3=cflat{3};
nth=zeros(size(mc3,1),size(mc3,2),2);
umx1c3=squeeze(Rinv(1,1)*mc3(:,:,4,:)+Rinv(1,2)*mc3(:,:,3,:)+Rinv(1,3)*mc3(:,:,2,:)+Rinv(1,4)*mc3(:,:,1,:));
umx2c3=squeeze(Rinv(2,1)*mc3(:,:,4,:)+Rinv(2,2)*mc3(:,:,3,:)+Rinv(2,3)*mc3(:,:,2,:)+Rinv(2,4)*mc3(:,:,1,:));
umx3c3=squeeze(Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:)+Rinv(3,4)*mc3(:,:,1,:));


col1=mean(umx1c3,3); col1(col1<0) = 0;
col2=mean(umx2c3,3); col2(col2<0) = 0; 
col3=mean(umx3c3,3); col3(col3<0) = 0;
lincoeff = 1;
settings.mc.lincoeff3 = lincoeff;
r3=(col3).^lincoeff;
g3=(col2).^lincoeff;
b3=(col1).^lincoeff;
r3=r3/max(r3(:));
% imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g3=g3/max(g3(:));
% imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b3=b3/max(b3(:));
% imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

merged3=cat(3,r3,g3,b3);
merged3=1.5*merged3;
merged3(merged3>1)=1;

% figure(6), 
% subplot(131);imshow(b3,[0 1]);colormap(myBmap);freezeColors;
% subplot(132);imshow(g3,[0 1]);colormap(myGmap);freezeColors;
% subplot(133);imshow(r3,[0 1]);colormap(myRmap);freezeColors;

figure(7);
image(merged3);
set(gca,'xtick',[],'ytick',[])
%% Linearization
disp('Linearization')
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
im1 = []; im2 = [] ; im3 = [];

for n=1:size(umx1c2,3)
%     disp(['Linearization un:',num2str(n)]);
    
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
% figure('units','normalized','outerposition',[0 0 1 1]); 
% subplot(131);imshow(mean(im1,3),[]);colormap(myBmap);freezeColors;title('channel 1')
% subplot(132);imshow(mean(im2,3),[]);colormap(myGmap);freezeColors;title('virtual channel') 
% subplot(133);imshow(mean(im3,3),[]);colormap(myRmap);freezeColors;title('channel 2')


%% Merging linearized deconvolved unmixed cumulants into final RGB image
bcgsub = 0.5;
r_lin=mean(im3,3).^lincoeff;
g_lin=mean(im2,3).^lincoeff;
b_lin=mean(im1,3).^lincoeff;
r_lin=r_lin/max(r_lin(:));
g_lin=g_lin/max(g_lin(:));
b_lin=b_lin/max(b_lin(:));
 

merged_lin=cat(3,r_lin,g_lin,b_lin);
% merged_lin=merged_lin;

figure(10);image(merged_lin)
figure(11);
subplot(131);imagesc(r_lin); title('Red')
subplot(132);imagesc(g_lin); title('Green')
subplot(133);imagesc(b_lin); title('Blue')

%%
% figure;
% image(merged);

figure(12)
subplot(131);
imshowpair(st1c,st2tc,'falsecolor');axis equal;axis tight;
title('STD overlay')
subplot(132);
imshow(mergedcum,[])
title('SOFI 2 unmixed')
subplot(133);
imshow(merged_lin,[]);
title('SOFI 2 lin dec unmixed')
set(gcf,'position',[465   555   1200   350])

figure(13)
subplot(131);
rf = zeros(size(r_lin,1),size(r_lin,2),3) ; rf(:,:,1) = r_lin;  rf(:,:,3) = r_lin; 
imshow(rf./max(rf(:)));
title('Red channel')
subplot(132);
gf = zeros(size(g_lin,1),size(g_lin,2),3) ; gf(:,:,2) = g_lin;  gf(:,:,1) = g_lin; 
imshow(gf./max(gf(:)));
title('Green channel')
subplot(133);
bf = zeros(size(b_lin,1),size(b_lin,2),3) ; bf(:,:,3) = b_lin;  bf(:,:,2) = b_lin; 
imshow(bf./max(bf(:)));
title('Blue channel')
set(gcf,'position',[465   555   1200   350])

%% saving 
disp('Saving results')
settings.io.outputpath = [pname,filesep,rFold];

% save std of stacks 
writeTIFF(st1c./max(st1c(:)),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_std_R'])
writeTIFF(st2tc./max(st2tc(:)),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_std_T'])
writeRGBTIFF(uint8(fake),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_coreg'])

 t =[];
% save SOFI cflat
% t(:,:,1) = mean(cflat{1}(:,:,1,:),4); t(:,:,2) = mean(cflat{1}(:,:,2,:),4);
% writeTIFF(t./max(t(:)),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_SOFI_cflat1']); t =[];
t(:,:,1) = mean(cflat{2}(:,:,1,:),4); t(:,:,2) = mean(cflat{2}(:,:,2,:),4); t(:,:,3) = mean(cflat{2}(:,:,3,:),4);
writeTIFF(t./max(t(:)),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_SOFI_cflat2']); t =[];
t(:,:,1) = mean(cflat{3}(:,:,1,:),4); t(:,:,2) = mean(cflat{3}(:,:,2,:),4); t(:,:,3) = mean(cflat{3}(:,:,3,:),4); t(:,:,4) = mean(cflat{3}(:,:,4,:),4);
writeTIFF(t./max(t(:)),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_SOFI_cflat3']); t =[];

% save unmixed 2nd
t(:,:,1) = r; t(:,:,2) = g; t(:,:,3) = b;
writeTIFF(t,[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_unmixed2_lin']); t =[];
writeRGBTIFF(uint8(255.*mergedcum),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_unmixed2_lin_overlay'])

% save unmixed 3rd
t(:,:,1) = r3; t(:,:,2) = g3; t(:,:,3) = b3;
writeTIFF(t,[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_unmixed3_lin']); t =[];
writeRGBTIFF(uint8(255.*merged3),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_unmixed3_lin_overlay'])

% save lin deconv
t(:,:,1) = r_lin; t(:,:,2) = g_lin; t(:,:,3) = b_lin;
writeTIFF(t,[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_unmixed2_dec']); t =[];
writeRGBTIFF(uint8(255.*merged_lin),[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_unmixed2_dec_overlay'])

saveFigure(12,[settings.io.pn,filesep,rFold],[tag,'_',rList(k).id,'_results'],'tif')
saveFigure(13,[settings.io.pn,filesep,rFold],[tag,'_',rList(k).id,'_results_RGB'],'tif')
saveSettingsTxt(settings,[settings.io.pn,filesep,rFold,filesep,tag,'_',rList(k).id,'_settings.txt'])

end




