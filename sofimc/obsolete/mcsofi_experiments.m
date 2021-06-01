clear all;

% impath = 'D:\Users\AzatSharipov\Measurements\20130201_C2C12_ABT20_Al647_formfixed\005\';
% impath = 'I:\asharipov\Measurements\20130304_C2C12_3colors\001_exp30ms\';
impath = 'G:\20120727_MultiColorMeasurements\Sample8_4colors\Sample8_4colors\Data17_RegBin\';
workingpath = 'D:\sgeissbu\Images\20130304_3colors\001\';

file1='c1.bin';
file2='c2.bin';

% load([impath 'settings.mat']);
%% load data

fid1=fopen([impath file1]);
fid2=fopen([impath file2]);
info=dir([impath file1]);

nx=512;%cam1.ROIPosition(3);
ny=512;%cam1.ROIPosition(4);

N=info.bytes/2/nx/ny;

k=0;
data1=uint16(zeros(ny,nx,N));
data2=data1;
imsum1=zeros(ny,nx);
imsum2=imsum1;
imxvar1=zeros(ny-1,nx);
imxvar2=imxvar1;
while k<N
    k=k+1
    data1(:,:,k)=fread(fid1,[ny,nx],'*uint16');
    data2(:,:,k)=fread(fid2,[ny,nx],'*uint16');
end
fclose(fid1);
fclose(fid2);

%%
subseqlength=500;
totlength=5000;
start=1;
Nss=floor((totlength-start+1)/subseqlength);
c2=zeros(2*size(data1,1)-1,2*size(data1,2)-1,3,Nss);
for ns=1:Nss
    ns
    fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
    d1=double(data1(:,:,fr));
    sd1=sum(d1,3);
    d2=double(data2(:,:,fr));
    sd2=sum(d2,3);
    
    disp('1. comb');
    c2(1:2:end,1:2:end,1,ns)=1/(subseqlength-1)*sum(d1(:,:,1:end-1).*d1(:,:,2:end),3)...
        -1/(subseqlength-1)^2*(sd1-d1(:,:,1)).*(sd1-d1(:,:,end));
    c2(1:2:end,1:2:end,2,ns)=1/(subseqlength-1)*sum(d1(:,:,1:end-1).*d2(:,:,2:end),3)...
        -1/(subseqlength-1)^2*(sd1-d1(:,:,1)).*(sd2-d2(:,:,end));
    c2(1:2:end,1:2:end,3,ns)=1/(subseqlength-1)*sum(d2(:,:,1:end-1).*d2(:,:,2:end),3)...
        -1/(subseqlength-1)^2*(sd2-d2(:,:,1)).*(sd2-d2(:,:,end));
    
    disp('2. comb');
    c2(2:2:end,1:2:end-1,1,ns)=1/(subseqlength-1)*sum(d1(2:end,1:end-1,1:end-1).*d1(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd1(2:end,1:end-1)-d1(2:end,1:end-1,1)).*(sd1(1:end-1,1:end-1)-d1(1:end-1,1:end-1,end));
    c2(2:2:end,1:2:end-1,2,ns)=1/(subseqlength-1)*sum(d1(2:end,1:end-1,1:end-1).*d2(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd1(2:end,1:end-1)-d1(2:end,1:end-1,1)).*(sd2(1:end-1,1:end-1)-d2(1:end-1,1:end-1,end));
    c2(2:2:end,1:2:end-1,3,ns)=1/(subseqlength-1)*sum(d2(2:end,1:end-1,1:end-1).*d2(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd2(2:end,1:end-1)-d2(2:end,1:end-1,1)).*(sd2(1:end-1,1:end-1)-d1(1:end-1,1:end-1,end));
    
    disp('3. comb');
    c2(2:2:end,2:2:end,1,ns)=1/(subseqlength-1)*sum(d1(2:end,2:end,1:end-1).*d1(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd1(2:end,2:end)-d1(2:end,2:end,1)).*(sd1(1:end-1,1:end-1)-d1(1:end-1,1:end-1,end));
    c2(2:2:end,2:2:end,2,ns)=1/(subseqlength-1)*sum(d1(2:end,2:end,1:end-1).*d2(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd1(2:end,2:end)-d1(2:end,2:end,1)).*(sd2(1:end-1,1:end-1)-d2(1:end-1,1:end-1,end));
    c2(2:2:end,2:2:end,3,ns)=1/(subseqlength-1)*sum(d2(2:end,2:end,1:end-1).*d2(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd2(2:end,2:end)-d2(2:end,2:end,1)).*(sd2(1:end-1,1:end-1)-d1(1:end-1,1:end-1,end));
    
    disp('4. comb');
    c2(1:2:end-1,2:2:end,1,ns)=1/(subseqlength-1)*sum(d1(1:end-1,2:end,1:end-1).*d1(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd1(1:end-1,2:end)-d1(1:end-1,2:end,1)).*(sd1(1:end-1,1:end-1)-d1(1:end-1,1:end-1,end));
    c2(1:2:end-1,2:2:end,2,ns)=1/(subseqlength-1)*sum(d1(1:end-1,2:end,1:end-1).*d2(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd1(1:end-1,2:end)-d1(1:end-1,2:end,1)).*(sd2(1:end-1,1:end-1)-d2(1:end-1,1:end-1,end));
    c2(1:2:end-1,2:2:end,3,ns)=1/(subseqlength-1)*sum(d2(1:end-1,2:end,1:end-1).*d2(1:end-1,1:end-1,2:end),3)...
        -1/(subseqlength-1)^2*(sd2(1:end-1,2:end)-d2(1:end-1,2:end,1)).*(sd2(1:end-1,1:end-1)-d1(1:end-1,1:end-1,end));
end

%% 2nd order histogram
U=abs([reshape(c2(100:900,100:900,1,2:end),[],1) reshape(c2(100:900,100:900,2,2:end),[],1)]);
V=abs([reshape(c2(100:900,100:900,2,2:end),[],1) reshape(c2(100:900,100:900,3,2:end),[],1)]);
nbins=500;
bins=linspace(0,pi/2,nbins);
N=histc(atan2(U(:,2),U(:,1)),bins);
M=histc(atan2(V(:,2),V(:,1)),bins);
x=tan(pi/2/(nbins-1)/2+bins(2:end-1));

maxn=max(N(1:end-1)+M(1:end-1));

% Nxy=hist3(U,[300 300])';
% Mxy=hist3(V,[300 300])';
% figure;imagesc(log10(Nxy(1:100,1:100)+Mxy(1:100,1:100)));axis xy;ylabel('TT');xlabel('RT');axis equal;axis tight;


figure;
semilogx(x,N(2:end-1));
hold on;
semilogx(x,M(2:end-1),'r');
semilogx(x,N(2:end-1)+M(2:end-1),'k');
grid on;
axis tight;

%% rgb camera
r=abs(mean(c2(100:900,100:900,3,2:end),4)).^(1/1);
g=abs(mean(c2(100:900,100:900,2,2:end),4)).^(1/1);
b=abs(mean(c2(100:900,100:900,1,2:end),4)).^(1/1);
r=0.7*r/max(r(:));
g=g/max(g(:));
b=1.7*b/max(b(:));

nth=zeros(size(r,1),size(r,2),2);
imwrite(uint8(2^8*cat(3,r/0.7,nth)),'sc2_3.png');
imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'sc2_2.png');
imwrite(uint8(2^8*cat(3,nth,b/1.3)),'sc2_1.png');

merged=cat(3,r,g,b);
merged=2*merged;
merged(merged>1)=1;
figure;
image(merged);
imwrite(uint8(2^8*merged),'merged_mixed.png');

%%
m1=mean(double(data1(50:450,50:450,1:100)),3);
m2=mean(double(data2(50:450,50:450,1:100)),3);
m1=m1-min(m1(:));
m2=m2-min(m2(:));
m1=m1/max(m1(:));
m2=m2/max(m2(:));
nth=zeros(size(m1,1),size(m1,2),2);
imwrite(uint8(2^8*cat(3,nth,m1)),'m1.png');
imwrite(uint8(2^8*cat(3,m2,nth)),'m2.png');
imwrite(uint8(2^8*cat(3,m2,zeros(size(m1)),m1)),'m12merged.png');
%% unmix 3 colors
R1=0.01;
R2=0.57;
R3=100;
R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1];
Rinv=R^(-1);

mc2=c2(100:900,100:900,:,1:end);
nth=zeros(size(mc2,1),size(mc2,2),2);
umx1c2=squeeze(Rinv(1,1)*mc2(:,:,3,:)+Rinv(1,2)*mc2(:,:,2,:)+Rinv(1,3)*mc2(:,:,1,:));
umx2c2=squeeze(Rinv(2,1)*mc2(:,:,3,:)+Rinv(2,2)*mc2(:,:,2,:)+Rinv(2,3)*mc2(:,:,1,:));
umx3c2=squeeze(Rinv(3,1)*mc2(:,:,3,:)+Rinv(3,2)*mc2(:,:,2,:)+Rinv(3,3)*mc2(:,:,1,:));

% figure;imagesc(umx1c2,[0 max(umx1c2(:))]);
% figure;imagesc(umx2c2,[0 max(umx2c2(:))]);
% figure;imagesc(umx3c2,[0 max(umx3c2(:))]);

col1=mean(umx1c2(:,:,2:end),3);
col2=mean(umx2c2(:,:,2:end),3);
col3=mean(umx3c2(:,:,2:end),3);
r=abs(col3).^(1/1);
g=abs(col2).^(1/1);
b=abs(col1).^(1/1);
r=0.7*r/max(r(:));
imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=1.2*b/max(b(:));
imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

merged=cat(3,r,g,b);
merged=2*merged;
merged(merged>1)=1;
figure;
image(merged);
imwrite(uint8(2^8*merged),'merged_unmixed.png');

%% distance factor determination color 1
U=abs([reshape(umx1c2(1:2:end-1,1:2:end-1,2:end),[],1) reshape(umx1c2(2:2:end,1:2:end-1,2:end),[],1)]);
V=abs([reshape(umx1c2(1:2:end-1,1:2:end-1,2:end),[],1) reshape(umx1c2(1:2:end-1,2:2:end,2:end),[],1)]);
W=abs([reshape(umx1c2(1:2:end-1,1:2:end-1,2:end),[],1) reshape(umx1c2(2:2:end,2:2:end,2:end),[],1)]);

% tmp=U(:,1);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% U(:,1)=tmp;
% tmp=U(:,2);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% U(:,2)=tmp;
% tmp=V(:,1);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% V(:,1)=tmp;
% tmp=V(:,2);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% V(:,2)=tmp;
% tmp=W(:,1);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% W(:,1)=tmp;
% tmp=W(:,2);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% W(:,2)=tmp;

nbins=100;
bins=linspace(0,pi/2,nbins);
N=histc(atan2(U(:,2),U(:,1)),bins);
M=histc(atan2(V(:,2),V(:,1)),bins);
O=histc(atan2(W(:,2),W(:,1)),bins);
x=tan(pi/2/(nbins-1)/2+bins(2:end-1));

figure;
semilogx(x,N(2:end-1));
hold on;
semilogx(x,M(2:end-1),'r');
semilogx(x,O(2:end-1),'g');
grid on;
axis tight;

% nbins=200;
% bins=linspace(0,pi/2,nbins);%atan(logspace(-2,2,200));
% 
% X=U;
% X(X<0)=0;
% tmp=X(:,1);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% X(:,1)=tmp;
% tmp=X(:,2);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% X(:,2)=tmp;
% 
% M=histc(atan2(X(:,2),X(:,1)),bins);
% 
% f=ezfit(bins(1:end-1),M(1:end-1),'gauss');
% df=tan(pi/2/(nbins-1)/2+f.m(3))
% 
% Mfit=f.m(1)*exp(-((f.x-f.m(3)).^2)./(2*f.m(2)^2));
% figure;plot(f.x,f.y);showfit(f);
% 
% X=V;
% X(X<0)=0;
% tmp=X(:,1);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% X(:,1)=tmp;
% tmp=X(:,2);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% X(:,2)=tmp;
% 
% M=histc(atan2(X(:,2),X(:,1)),bins);
% 
% f=ezfit(bins(1:end-1),M(1:end-1),'gauss');
% df=tan(pi/2/(nbins-1)/2+f.m(3))
% 
% Mfit=f.m(1)*exp(-((f.x-f.m(3)).^2)./(2*f.m(2)^2));
% figure;plot(f.x,f.y);showfit(f);
% 
% X=W;
% X(X<0)=0;
% tmp=X(:,1);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% X(:,1)=tmp;
% tmp=X(:,2);
% tmp(tmp<0.02*max(tmp(:)))=NaN;
% X(:,2)=tmp;
% 
% M=histc(atan2(X(:,2),X(:,1)),bins);
% 
% f=ezfit(bins(1:end-1),M(1:end-1),'gauss');
% df=tan(pi/2/(nbins-1)/2+f.m(3))
% 
% Mfit=f.m(1)*exp(-((f.x-f.m(3)).^2)./(2*f.m(2)^2));
% figure;plot(f.x,f.y);showfit(f);

%% flatten color 1
col1=mean(umx1c2(:,:,2:end),3);
figure;imagesc(col1,[0 max(col1(:))]);colormap(fireprint(256));axis image;

d=1.1;
col1(2:2:end,1:2:end-1)=col1(2:2:end,1:2:end-1)/d;
umx1c2(2:2:end,1:2:end-1,:)=umx1c2(2:2:end,1:2:end-1,:)/d;
d=1.1;
col1(1:2:end-1,2:2:end)=col1(1:2:end-1,2:2:end)/d;
umx1c2(1:2:end-1,2:2:end,:)=umx1c2(1:2:end-1,2:2:end,:)/d;
d=1.25;
col1(2:2:end,2:2:end)=col1(2:2:end,2:2:end)/d;
umx1c2(2:2:end,2:2:end,:)=umx1c2(2:2:end,2:2:end,:)/d;
figure;imagesc(col1,[0 max(col1(:))]);colormap(fireprint(256));axis image;

%% distance factor determination color 2
U=abs([reshape(umx2c2(1:2:end-1,1:2:end-1,2:end),[],1) reshape(umx2c2(2:2:end,1:2:end-1,2:end),[],1)]);
V=abs([reshape(umx2c2(1:2:end-1,1:2:end-1,2:end),[],1) reshape(umx2c2(1:2:end-1,2:2:end,2:end),[],1)]);
W=abs([reshape(umx2c2(1:2:end-1,1:2:end-1,2:end),[],1) reshape(umx2c2(2:2:end,2:2:end,2:end),[],1)]);

tmp=U(:,1);
tmp(tmp<0.02*max(tmp(:)))=NaN;
U(:,1)=tmp;
tmp=U(:,2);
tmp(tmp<0.02*max(tmp(:)))=NaN;
U(:,2)=tmp;
tmp=V(:,1);
tmp(tmp<0.02*max(tmp(:)))=NaN;
V(:,1)=tmp;
tmp=V(:,2);
tmp(tmp<0.02*max(tmp(:)))=NaN;
V(:,2)=tmp;
tmp=W(:,1);
tmp(tmp<0.02*max(tmp(:)))=NaN;
W(:,1)=tmp;
tmp=W(:,2);
tmp(tmp<0.02*max(tmp(:)))=NaN;
W(:,2)=tmp;

nbins=100;
bins=linspace(0,pi/2,nbins);
N=histc(atan2(U(:,2),U(:,1)),bins);
M=histc(atan2(V(:,2),V(:,1)),bins);
O=histc(atan2(W(:,2),W(:,1)),bins);
x=tan(pi/2/(nbins-1)/2+bins(2:end-1));

figure;
semilogx(x,N(2:end-1));
hold on;
semilogx(x,M(2:end-1),'r');
semilogx(x,O(2:end-1),'g');
grid on;
axis tight;

%% flatten color 2
col2=mean(umx2c2(:,:,2:end),3);
figure;imagesc(col2,[0 max(col2(:))]);colormap(fireprint(256));axis image;

d=1.1;
col2(2:2:end,1:2:end-1)=col2(2:2:end,1:2:end-1)/d;
umx2c2(2:2:end,1:2:end-1,:)=umx2c2(2:2:end,1:2:end-1,:)/d;
d=1.05;
col2(1:2:end-1,2:2:end)=col2(1:2:end-1,2:2:end)/d;
umx2c2(1:2:end-1,2:2:end,:)=umx2c2(1:2:end-1,2:2:end,:)/d;
d=1.2;
col2(2:2:end,2:2:end)=col2(2:2:end,2:2:end)/d;
umx2c2(2:2:end,2:2:end,:)=umx2c2(2:2:end,2:2:end,:)/d;
figure;imagesc(col2,[0 max(col2(:))]);colormap(fireprint(256));axis image;

%% distance factor determination color 3
U=abs([reshape(umx3c2(1:2:end-1,1:2:end-1,3:end),[],1) reshape(umx3c2(2:2:end,1:2:end-1,3:end),[],1)]);
V=abs([reshape(umx3c2(1:2:end-1,1:2:end-1,3:end),[],1) reshape(umx3c2(1:2:end-1,2:2:end,3:end),[],1)]);
W=abs([reshape(umx3c2(1:2:end-1,1:2:end-1,3:end),[],1) reshape(umx3c2(2:2:end,2:2:end,3:end),[],1)]);

tmp=U(:,1);
tmp(tmp<0.02*max(tmp(:)))=NaN;
U(:,1)=tmp;
tmp=U(:,2);
tmp(tmp<0.02*max(tmp(:)))=NaN;
U(:,2)=tmp;
tmp=V(:,1);
tmp(tmp<0.02*max(tmp(:)))=NaN;
V(:,1)=tmp;
tmp=V(:,2);
tmp(tmp<0.02*max(tmp(:)))=NaN;
V(:,2)=tmp;
tmp=W(:,1);
tmp(tmp<0.02*max(tmp(:)))=NaN;
W(:,1)=tmp;
tmp=W(:,2);
tmp(tmp<0.02*max(tmp(:)))=NaN;
W(:,2)=tmp;

nbins=50;
bins=linspace(0,pi/2,nbins);
N=histc(atan2(U(:,2),U(:,1)),bins);
M=histc(atan2(V(:,2),V(:,1)),bins);
O=histc(atan2(W(:,2),W(:,1)),bins);
x=tan(pi/2/(nbins-1)/2+bins(2:end-1));

figure;
semilogx(x,N(2:end-1));
hold on;
semilogx(x,M(2:end-1),'r');
semilogx(x,O(2:end-1),'g');
grid on;
axis tight;

%% flatten color 3
col3=mean(umx3c2(:,:,3:end),3);
figure;imagesc(col3,[0 max(col3(:))]);colormap(fireprint(256));axis image;

d=1.12;
col3(2:2:end,1:2:end-1)=col3(2:2:end,1:2:end-1)/d;
umx3c2(2:2:end,1:2:end-1,:)=umx3c2(2:2:end,1:2:end-1,:)/d;
d=1.12;
col3(1:2:end-1,2:2:end)=col3(1:2:end-1,2:2:end)/d;
umx3c2(1:2:end-1,2:2:end,:)=umx3c2(1:2:end-1,2:2:end,:)/d;
d=1.13;
col3(2:2:end,2:2:end)=col3(2:2:end,2:2:end)/d;
umx3c2(2:2:end,2:2:end,:)=umx3c2(2:2:end,2:2:end,:)/d;
figure;imagesc(col3,[0 max(col3(:))]);colormap(fireprint(256));axis image;

%% save flattened colors
col1(col1<0)=0;
col2(col2<0)=0;
col3(col3<0)=0;
imwrite(uint8(2^8*abs(col1)/max(abs(col1(:)))),'col1.png');
imwrite(uint8(2^8*abs(col2)/max(abs(col2(:)))),'col2.png');
imwrite(uint8(2^8*abs(col3)/max(abs(col3(:)))),'col3.png');

r=abs(col3).^(1/1);
g=abs(col2).^(1/1);
b=abs(col1).^(1/1);
r=0.7*r/max(r(:));
imwrite(uint8(2^8*cat(3,r/0.7,nth)),'col3.png');
g=g/max(g(:));
imwrite(uint8(2^8*cat(3,zeros(size(g)),g,zeros(size(g)))),'col2.png');
b=1.2*b/max(b(:));
imwrite(uint8(2^8*cat(3,nth,b/1.2)),'col1.png');

merged=cat(3,r,g,b);
merged=2*merged;
merged(merged>1)=1;
figure;
image(merged);
imwrite(uint8(2^8*merged),'merged_unmixed.png');

%% deconvolve
fwhm3=4.2;
fwhm2=4;
fwhm1=3.8;

dccol3=umx3c2;
dccol3(dccol3<0.002*max(dccol3(:)))=0;
for ns=1:size(dccol3,3)
    ns
    scale=max(max(dccol3(:,:,ns)));
    dccol3(:,:,ns)=scale*deconvlucy(dccol3(:,:,ns)/scale,fspecial('Gaussian',15,fwhm3/sqrt(8*log(2))),10);
end
figure;imagesc(abs(mean(dccol3(:,:,2:end),3)));colormap(fireprint(256));axis image;

dccol2=umx2c2;
dccol2(dccol2<0)=0;
for ns=1:size(dccol2,3)
    ns
    scale=max(max(dccol2(:,:,ns)));
    dccol2(:,:,ns)=scale*deconvlucy(dccol2(:,:,ns)/scale,fspecial('Gaussian',15,fwhm2/sqrt(8*log(2))),10);
end
figure;imagesc(abs(mean(dccol2(:,:,2:end),3)));colormap(fireprint(256));axis image;

dccol1=umx1c2;
dccol1(dccol1<0)=0;
for ns=1:size(dccol1,3)
    ns
    scale=max(max(dccol1(:,:,ns)));
    dccol1(:,:,ns)=scale*deconvlucy(dccol1(:,:,ns)/scale,fspecial('Gaussian',15,fwhm1/sqrt(8*log(2))),10);
end
figure;imagesc(abs(mean(dccol1(:,:,2:end),3)));colormap(fireprint(256));axis image;
%%
umx1c2=abs(umx1c2)/max(abs(umx1c2(:)));
umx2c2=abs(umx2c2)/max(abs(umx2c2(:)));
umx3c2=abs(umx3c2)/max(abs(umx3c2(:)));
col=isolum_cbs(3);
col=flipud(col);
tmp=0;
figure;
cim=(repmat(umx1c2,[1 1 3]).*repmat(reshape(col(1,:),[1 1 3]),[size(umx1c2,1) size(umx1c2,2) 1]));
image(cim/max(cim(:)));
axis image;
axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['col1.png']);
tmp=tmp+1/3*cim/max(cim(:));
figure;
cim=(repmat(umx2c2,[1 1 3]).*repmat(reshape(col(2,:),[1 1 3]),[size(umx1c2,1) size(umx1c2,2) 1]));
image(cim/max(cim(:)));
axis image;
axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['col2.png']);
tmp=tmp+1/3*cim/max(cim(:));
figure;
cim=(repmat(umx3c2,[1 1 3]).*repmat(reshape(col(3,:),[1 1 3]),[size(umx1c2,1) size(umx1c2,2) 1]));
image(cim/max(cim(:)));
axis image;
axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['col3.png']);
tmp=tmp+1/3*cim/max(cim(:));

tmp=tmp/max(tmp(:));
figure;
image(tmp);
imwrite(uint8(2^8*tmp),['mergedcolorsc1.png']);
axis image;
axis off;

%%
subseqlength=1000;
totlength=5000;
start=1;
pixelroi=[];
orders=1:4;
fwhmest=4;%[]%2*sqrt(8*log(2));

Nss=floor((totlength-start+1)/subseqlength);

c=cell(numel(orders),1);
for m=orders
    c{m}=[];
end
bc=c;
fwhm=zeros(Nss,1);
for n=1:Nss
    disp(n);
    [sofi,grids]=sofiCumulants2D(data2,start+(n-1)*subseqlength,subseqlength,pixelroi,orders);
%     [sofi,fwhm(n)]=sofiFlatten(fwhmest,sofi,grids);
    for m=orders
        c{m}(:,:,end+1)=sofi{m};
    end
end
for m=orders
    c{m}(:,:,1)=[];
end

%%
% extracting the center of gravities of individual fluorophores
sys=struct;
out=struct;
out.ru=5;
sys.bg=1;

bgth=-1;
logsize=2.1;
alol=15;
aupl=25;

cog_fix=[];
cog_mov=[];
for n=2:size(c{1},3)
    n
    im_fix=abs(cdata2{3}(1:3:end,1:3:end,n)).^(1/3);
    im_mov=abs(cdata1{3}(1:3:end,1:3:end,n)).^(1/3);

    im_fix=im_fix.*mask;
    im_mov=im_mov.*mask;
%     figure;imagesc(im_fix);
    out=hriSegmentation(double(im_fix),bgth,logsize,out);
    bgmfix=out.bgmap;
%     figure;imagesc(out.bgmap);
    out=hriFilterSegments(double(im_fix),aupl,alol,sys,out);
    cog_fix=cat(1,cog_fix,[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate]);
    out=hriSegmentation(double(im_mov),bgth,logsize,out);
    bgmmov=out.bgmap;
    out=hriFilterSegments(double(im_mov),aupl,alol,sys,out);
    cog_mov=cat(1,cog_mov,[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate]);
end

% identify corresponding beads in the two channels
% spatial cross-correlation algorithm to determine shift of coordinates
res=normxcorr2(bgmfix,bgmmov);    
[~,m_corr]=max(res(:));
[mx,my]=ind2sub(size(res),m_corr);
mx=-mx+size(im_fix,1)
my=-my+size(im_fix,2)

cog_common_fix=zeros(1,2);
cog_common_mov=zeros(1,2);

%identify corresponding center of gravities
pixel_tolerance=3;
for k=1:size(cog_fix,1)
    d=(cog_mov-repmat(cog_fix(k,:)-[mx my],[size(cog_mov,1) 1]));
    d=sqrt(d(:,1).^2+d(:,2).^2);
    [p,ind]=min(d);
    if(p<pixel_tolerance)
        cog_common_fix(end+1,:)=cog_fix(k,:);
        cog_common_mov(end+1,:)=cog_mov(ind,:);
    end
end

cog_common_fix(1,:)=[];  %first element is irrelevant (0,0)
cog_common_mov(1,:)=[];

figure(20);
imshowpair(double(im_fix), double(im_mov),'Scaling','joint'); %green: first image, magenta: second image
hold on;
scatter(cog_fix(:,2),cog_fix(:,1),'v','MarkerEdgeColor',[0 1 0]);
scatter(cog_mov(:,2),cog_mov(:,1),'v','MarkerEdgeColor',[1 0 1]);
scatter(cog_common_fix(:,2),cog_common_fix(:,1),'o','MarkerEdgeColor',[0 1 0]);
scatter(cog_common_mov(:,2),cog_common_mov(:,1),'o','MarkerEdgeColor',[1 0 1]);

d=(cog_common_mov+repmat([mx my],[size(cog_common_mov,1) 1]))-cog_common_fix;
disp(['avg error of COG coordinates using pure displacement at pixel level: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);

% affine transformation
ccm=fliplr(cog_common_mov);
ccf=fliplr(cog_common_fix);
tf{m+1}=cp2tform(ccm,ccf,'affine');

[x,y] = tformfwd(tf{m+1},cog_common_mov(:,2),cog_common_mov(:,1));

figure(21);
clf(21);
scatter(cog_common_fix(:,2),cog_common_fix(:,1),'o','MarkerEdgeColor',[0 1 0]);
hold on;
scatter(x,y,'o','MarkerEdgeColor',[1 0 1]);

d=([y x]-cog_common_fix);
disp(['avg error of COG coordinates using affine transformation: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);

moving_reg = imtransform(im_mov, tf{m+1}, 'XData', [1 size(im_fix,2)], 'YData', [1 size(im_fix,1)], ...
   'Size', size(im_fix));

% [x,y] = tformfwd(tfexpcp{m+1},cog_common_mov(:,2),cog_common_mov(:,1));
% d=([y x]-cog_common_fix);
% disp(['avg error of COG coordinates using affine transformation from later calibration: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);

figure(22);
imshowpair(log(double(im_fix)+100), log(double(moving_reg)+100),'Scaling','joint'); %green: first image, magenta: second image

% T=1;
% for mi=0:m
%     T=T*tf{mi+1}.tdata.T;
% end
% tftmp=maketform('affine',T);
% 
% transformed_stacks(:,:,:,m+2) = imtransform(img_seq_ch_moving, tftmp, 'XData', [1 size(im_fix,2)], 'YData', [1 size(im_fix,1)], ...
%    'Size', size(im_fix));
% figure;
% imshowpair(log(double(im0)+100), log(double(moving_reg)+100),'Scaling','joint'); %green: first image, magenta: second image

%% Calculate raw cross-cumulants
subseqlength=1000;
totlength=5000;
start=1;

Nss=floor((totlength-start+1)/subseqlength);

orders=1:4;
tf=tform;

% if ~exist([workingpath 'evalss' num2str(subseqlength)])
%     mkdir([workingpath 'evalss' num2str(subseqlength)]);
% end
% evalpath=[workingpath 'evalss' num2str(subseqlength) '\'];
c=cell(4,1);
c2top=[];
for ns=1:Nss
    ns
    
    fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
    im_fix=double(data1(:,:,fr));
    im_mov=double(data2(:,:,fr));


    tic;
%         mov_reg = permute(imtransform(permute(im_mov,[2 1 3]), tf{m+1}, 'XData', [1 size(im_fix,2)], 'YData', [1 size(im_fix,1)], ...
%             'Size', size(im_fix)),[2 1 3]);
%         mov_reg = imtransform(im_mov, tfmavg{m+1}, 'XData', [1 size(im_fix,2)], 'YData', [1 size(im_fix,1)], ...
%             'Size', size(im_fix));
    T=tf.tdata.T;
    mov_reg = imtransform(im_mov,maketform('affine',T), 'XData', [1 size(im_fix,2)], 'YData', [1 size(im_fix,1)], ...
        'Size', size(im_fix));
    toc

    [sofi,grids]=mcSofiCumulants(permute(cat(4,im_fix,mov_reg),[1 2 4 3]));
    for io=orders
        c{io}=cat(3,c{io},sofi{io}(:,:,1:end-1));
    end
    [sofi,grid2d]=mcSofiCumulantsSingle(im_mov);
    for io=orders
        c{io}=cat(3,c{io},sofi{io});
    end
end
for io=orders
    c{io}=reshape(c{io},size(c{io},1),size(c{io},2),[],Nss);
end

%%
%% histograms
U=abs([reshape(c{3}(:,:,1,:),[],1) reshape(c{3}(:,:,2,:),[],1)]);
V=abs([reshape(c{3}(:,:,1,:),[],1) reshape(c{3}(:,:,3,:),[],1)]);
W=abs([reshape(c{3}(:,:,1,:),[],1) reshape(c{3}(:,:,4,:),[],1)]);
X=abs([reshape(c{3}(:,:,2,:),[],1) reshape(c{3}(:,:,3,:),[],1)]);
Y=abs([reshape(c{3}(:,:,2,:),[],1) reshape(c{3}(:,:,4,:),[],1)]);
Z=abs([reshape(c{3}(:,:,3,:),[],1) reshape(c{3}(:,:,4,:),[],1)]);

% U=abs([reshape(c{3}(1:3:end,1:3:end,1,:),[],1) reshape(c{3}(1:3:end,1:3:end,2,:),[],1)]);
% V=abs([reshape(c{3}(1:3:end,1:3:end,1,:),[],1) reshape(c{3}(1:3:end,1:3:end,3,:),[],1)]);
% W=abs([reshape(c{3}(1:3:end,1:3:end,1,:),[],1) reshape(c{3}(1:3:end,1:3:end,4,:),[],1)]);
% X=abs([reshape(c{3}(1:3:end,1:3:end,2,:),[],1) reshape(c{3}(1:3:end,1:3:end,3,:),[],1)]);
% Y=abs([reshape(c{3}(1:3:end,1:3:end,2,:),[],1) reshape(c{3}(1:3:end,1:3:end,4,:),[],1)]);
% Z=abs([reshape(c{3}(1:3:end,1:3:end,3,:),[],1) reshape(c{3}(1:3:end,1:3:end,4,:),[],1)]);

% Nxy=hist3(U,[300 300])';
% Mxy=hist3(X,[300 300])';
% Oxy=hist3(Z,[300 300])';
% figure;imagesc(log10(Mxy(1:100,1:100)));axis xy;ylabel('RTT');xlabel('RRT');axis equal;axis tight;

nbins=200;
bins=linspace(0,pi/2,nbins);%atan(logspace(-2,2,200));
N=histc(atan2(U(:,2),U(:,1)),bins);
M=histc(atan2(X(:,2),X(:,1)),bins);
O=histc(atan2(Z(:,2),Z(:,1)),bins);
P=histc(atan2(V(:,2),V(:,1)),bins);
Q=histc(atan2(W(:,2),W(:,1)),bins);
% [N,binn]=hist(atan2(U(:,2),U(:,1)),200);
% [M,binm]=hist(atan2(X(:,2),X(:,1)),200);
% [O,bino]=hist(atan2(Z(:,2),Z(:,1)),200);
% [P,binp]=hist(atan2(V(:,2),V(:,1)),200);
% [Q,binq]=hist(atan2(W(:,2),W(:,1)),200);

x=tan(pi/2/(nbins-1)/2+bins(1:end-1));
maxn=max(N(1:end-1)+M(1:end-1)+O(1:end-1));

figure(4);
clf(4);
semilogx(x,N(1:end-1),'LineWidth',2);
hold on;
semilogx(x,M(1:end-1),'r','LineWidth',2);
semilogx(x,O(1:end-1),'g','LineWidth',2);
% semilogx((tan(binp)).^(1/2),P,'k');
% semilogx((tan(binq)).^(1/3),Q,'m');
semilogx(x,N(1:end-1)+M(1:end-1)+O(1:end-1),'k','LineWidth',2);
legend('RRT/RRR','RTT/RRT','TTT/RTT','RRT/RRR+RTT/RRT+TTT/RTT','R_1 Alexa 488','R_2 Rhodamine 6G','R_3 Alexa 568','R_4 Alexa 647');
grid on;
axis tight;

%% flatten cumulants
ctmp=c;
tmp=cell(numel(orders),1);
for ns=1:size(c{1},4)
    ns
    for io=orders
        tmp{io}=ctmp{io}(:,:,:,ns);
    end
    [tmp,fit]=sofiFlatten(4,tmp,grid2d);
%     fit
    for io=orders
        ctmp{io}(:,:,:,ns)=tmp{io};
    end
end