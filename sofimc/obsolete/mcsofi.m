% path='D:\sgeissbu\Images\MultiColorSim_4colors\';
path='C:\Users\Lukestom\Documents\ÈVUT\Ph.D. III.semestr\SOFI method\Original\MulticolorSOFI\MultiColorSim_4colors';
% clear all;
%%
subseqlength=1000;
totlength=10000;
start=1;
pixelroi=[];
orders=1:4;

Nss=floor((totlength-start+1)/subseqlength);

c=cell(numel(orders),1);
for m=orders
    c{m}=[];
end
bc=c;
fwhm=zeros(Nss,1);
for n=1:Nss
    disp(n);
    [sofi,grid]=mcSofiCumulants(stack,start+(n-1)*subseqlength,subseqlength,pixelroi,orders);
%     [sofi,fwhm(n)]=sofiFlatten(fwhmest,sofi,grid);
    for m=orders
        c{m}(:,:,:,end+1)=sofi{m};
    end
end
for m=orders
    c{m}(:,:,:,1)=[];
end

%% eval sum stack

subseqlength=1000;
totlength=10000;
start=1;
pixelroi=[];
orders=1:4;
fwhmest=2*sqrt(8*log(2));%[]%2*sqrt(8*log(2));

Nss=floor((totlength-start+1)/subseqlength);

c=cell(numel(orders),1);

bc=c;
fwhm=zeros(Nss,1);
for n=1:Nss
    disp(n);
    [sofi,grid]=sofiCumulants2D(sumstack,start+(n-1)*subseqlength,subseqlength,pixelroi,orders);
    [sofi,fwhm(n)]=sofiFlatten(fwhmest,sofi,grid);
    sofi=sofiLinearize(sofi,fwhmest,0.01,1:4,100);
    for m=orders
        c{m}=cat(3,c{m},sofi{m});
    end
end

%%
fwhm=2*sqrt(8*log(2));
R=[0.10,0.51,6.03,43.1];
Rm=[R(1)^3 R(2)^3 R(3)^3 R(4)^3; R(1)^2 R(2)^2 R(3)^2 R(4)^2; R(1) R(2) R(3) R(4); 1 1 1 1];
Rinv=Rm^(-1);
colors=cell(numel(orders),1);
colorscflat=colors;
colorscraw=colors;
grid2d=mcSofiGrids(3,2);
for n=1:Nss
    n
    sofi=cell(numel(orders),1);
    for m=orders
        sofi{m}=c{m}(:,:,:,n);
    end
    %unmixing 4 colors using 3rd order
    umx1c3=Rinv(1,1)*sofi{3}(:,:,4)+Rinv(1,2)*sofi{3}(:,:,3)+Rinv(1,3)*sofi{3}(:,:,2)+Rinv(1,4)*sofi{3}(:,:,1);
    umx2c3=Rinv(2,1)*sofi{3}(:,:,4)+Rinv(2,2)*sofi{3}(:,:,3)+Rinv(2,3)*sofi{3}(:,:,2)+Rinv(2,4)*sofi{3}(:,:,1);
    umx3c3=Rinv(3,1)*sofi{3}(:,:,4)+Rinv(3,2)*sofi{3}(:,:,3)+Rinv(3,3)*sofi{3}(:,:,2)+Rinv(3,4)*sofi{3}(:,:,1);
    umx4c3=Rinv(4,1)*sofi{3}(:,:,4)+Rinv(4,2)*sofi{3}(:,:,3)+Rinv(4,3)*sofi{3}(:,:,2)+Rinv(4,4)*sofi{3}(:,:,1);

    sofi=cell(3,1);
    sofi{3}=umx1c3;
    colorscraw{1}=cat(3,colorscraw{1},sofi{3});
    sofi=sofiFlatten(fwhm*488/647,sofi,grid2d,3);
    colorscflat{1}=cat(3,colorscflat{1},sofi{3});
    sofi=sofiLinearize(sofi,fwhm*488/647,0.01,3,100);
    psf=fspecial('Gaussian',13,1*488/647);
    sofi{3}(sofi{3}<0.2*max(sofi{3}(:)))=0;
    sofi{3}=conv2(sofi{3},psf,'same');
    umx1c3=sofi{3};
    sofi{3}=umx2c3;
    colorscraw{2}=cat(3,colorscraw{2},sofi{3});
    sofi=sofiFlatten(fwhm*520/647,sofi,grid2d,3);
    colorscflat{2}=cat(3,colorscflat{2},sofi{3});
    sofi=sofiLinearize(sofi,fwhm*520/647,0.01,3,100);
    psf=fspecial('Gaussian',13,1*520/647);
    sofi{3}(sofi{3}<0.2*max(sofi{3}(:)))=0;
    sofi{3}=conv2(sofi{3},psf,'same');
    umx2c3=sofi{3};
    sofi{3}=umx3c3;
    colorscraw{3}=cat(3,colorscraw{3},sofi{3});
    sofi=sofiFlatten(fwhm*568/647,sofi,grid2d,3);
    colorscflat{3}=cat(3,colorscflat{3},sofi{3});
    sofi=sofiLinearize(sofi,fwhm*568/647,0.01,3,100);
    psf=fspecial('Gaussian',13,1*568/647);
    sofi{3}(sofi{3}<0.2*max(sofi{3}(:)))=0;
    sofi{3}=conv2(sofi{3},psf,'same');
    umx3c3=sofi{3};
    sofi{3}=umx4c3;
    colorscraw{4}=cat(3,colorscraw{4},sofi{3});
    sofi=sofiFlatten(fwhm,sofi,grid2d,3);
    colorscflat{4}=cat(3,colorscflat{4},sofi{3});
    sofi=sofiLinearize(sofi,fwhm,0.01,3,100);
    psf=fspecial('Gaussian',13,1);
    sofi{3}(sofi{3}<0.2*max(sofi{3}(:)))=0;
    sofi{3}=conv2(sofi{3},psf,'same');
    umx4c3=sofi{3};
    
    colors{1}=cat(3,colors{1},umx1c3);
    colors{2}=cat(3,colors{2},umx2c3);
    colors{3}=cat(3,colors{3},umx3c3);
    colors{4}=cat(3,colors{4},umx4c3);
end

%%
mcolors=colors;
for m=1:length(colors)
    mcolors{m}=mean(colors{m},3);
    mcolors{m}=mcolors{m}/max(mcolors{m}(:));
end
col=isolum_cbs(length(R));
col=flipud(col);
% alphamap=abs(out.bsofi);%abs((sofi{4}(:,:,1)+sofi{4}(:,:,end))/2);
% alphamap=alphamap/max(alphamap(:));
tmp=0;
for n=1:length(R)
    figure;
    cim=(repmat(mcolors{n},[1 1 3]).*repmat(reshape(col(n,:),[1 1 3]),[size(mcolors{1},1) size(mcolors{1},2) 1]));
    image(cim/max(cim(:)));
    axis image;
    axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['col' num2str(n) 'bc3' '.png']);
    
    tmp=tmp+1/length(R)*cim;
end
% tmp(:,:,1)=tmp(:,:,1)./mean(reshape(cell2mat(mcolors)',size(mcolors{1},1),size(mcolors{1},2),[]),3)';
% tmp(:,:,2)=tmp(:,:,2)./mean(reshape(cell2mat(mcolors)',size(mcolors{1},1),size(mcolors{1},2),[]),3)';
% tmp(:,:,3)=tmp(:,:,3)./mean(reshape(cell2mat(mcolors)',size(mcolors{1},1),size(mcolors{1},2),[]),3)';
tmp=1.1*tmp/max(tmp(:));
tmp(tmp>1)=1;
figure;
image(tmp);
imwrite(uint8(2^8*tmp),['mergedcolorsbc3.png']);
axis image;
axis off;

%%
col=isolum_cbs(length(R));
col=flipud(col);
mcolors=colorscflat;
for m=1:length(colorscflat)
    mcolors{m}=abs(mean(colorscflat{m},3));
    mcolors{m}=mcolors{m}/max(mcolors{m}(:));
end
tmp=0;
for n=1:length(R)
    figure;
    cim=(repmat(mcolors{n},[1 1 3]).*repmat(reshape(col(n,:),[1 1 3]),[size(mcolors{1},1) size(mcolors{1},2) 1]));
    image(cim/max(cim(:)));
    axis image;
    axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['col' num2str(n) 'c3' '.png']);
    tmp=tmp+1/length(R)*cim;
end
tmp=tmp/max(tmp(:));
figure;
image(tmp);
imwrite(uint8(2^8*tmp),['mergedcolorsc3.png']);
axis image;
axis off;

%%
col=isolum_cbs(length(R));
col=flipud(col);
mcolors=colorscraw;
for m=1:length(colorscraw)
    mcolors{m}=abs(mean(colorscraw{m},3));
    mcolors{m}=mcolors{m}/max(mcolors{m}(:));
end
tmp=0;
for n=1:length(R)
    figure;
    cim=(repmat(mcolors{n},[1 1 3]).*repmat(reshape(col(n,:),[1 1 3]),[size(mcolors{1},1) size(mcolors{1},2) 1]));
    image(cim/max(cim(:)));
    axis image;
    axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['col' num2str(n) 'craw3' '.png']);
    tmp=tmp+1/length(R)*cim;
end
tmp=tmp/max(tmp(:));
figure;
image(tmp);
imwrite(uint8(2^8*tmp),['mergedcolorscraw3.png']);
axis image;
axis off;

%%
col=isolum_cbs(size(c{1},3));
col=flipud(col);
tmp=0;
for n=1:size(c{1},3)
    mc1=mean(c{1}(:,:,n,:),4);
    figure;
    cim=(repmat(mc1,[1 1 3]).*repmat(reshape(col(n,:),[1 1 3]),[size(mc1,1) size(mc1,2) 1]));
    image(cim/max(cim(:)));
    axis image;
    axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['mch' num2str(n) '.png']);
    tmp=tmp+1/length(R)*cim;
end
tmp=tmp/max(tmp(:));
figure;
image(tmp);
imwrite(uint8(2^8*tmp),['mergedcolorsc1.png']);
axis image;
axis off;

%% colorcam c2
mc2=mean(c{2},4);
r=abs(mc2(:,:,1))/max(max(abs(mc2(:,:,1))));
g=abs(mc2(:,:,2))/max(max(abs(mc2(:,:,2))));
b=abs(mc2(:,:,3))/max(max(abs(mc2(:,:,3))));
figure;
image(cat(3,r(2:2:end,1:2:end),g(2:2:end,1:2:end),b(2:2:end,1:2:end)));
%%
col=isolum_cbs(size(c{3},3));
col=flipud(col);
% for n=1:size(c{3},3)
%     tmp=abs(mean(c{3}(:,:,n,:),4));
%     tmp=tmp/max(tmp(:));
%     imwrite(uint8(2^8*tmp),['c3' repmat('t',1,n-1) repmat('r',1,size(c{3},3)-n) '.png']);
%     figure;imagesc(mean(c{3}(:,:,n,:),4));axis image;colormap(gray);axis off;
% end

tmp=0;
for n=1:size(c{3},3)
    mc=abs(mean(c{3}(:,:,n,:),4));
    figure;
    cim=(repmat(mc,[1 1 3]).*repmat(reshape(col(n,:),[1 1 3]),[size(mc,1) size(mc,2) 1]));
    image(cim/max(cim(:)));
    axis image;
    axis off;
    imwrite(uint8(2^8*cim/max(cim(:))),['c3' repmat('t',1,n-1) repmat('r',1,size(c{3},3)-n) '.png']);
    tmp=tmp+1/size(c{3},3)*cim;
end
tmp=tmp/max(tmp(:));
figure;
image(tmp);
imwrite(uint8(2^8*tmp),['mergedxc3.png']);
axis image;
axis off;
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
line([R(1) R(1)],[0 maxn],'Linestyle','--','Color',col(1,:),'LineWidth',4);
line([R(2) R(2)],[0 maxn],'linestyle','--','Color',col(2,:),'LineWidth',4);
line([R(3) R(3)],[0 maxn],'linestyle','--','Color',col(3,:),'LineWidth',4);
line([R(4) R(4)],[0 maxn],'linestyle','--','Color',col(4,:),'LineWidth',4);
legend('RRT/RRR','RTT/RRT','TTT/RTT','RRT/RRR+RTT/RRT+TTT/RTT','R_1 Alexa 488','R_2 Rhodamine 6G','R_3 Alexa 568','R_4 Alexa 647');
grid on;
axis tight;