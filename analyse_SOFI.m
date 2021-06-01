%% Allows the use to load multiple sofi data and to perform various analysis on top

% pick any SOFI data
addpath(genpath('utils'))
rpath = 'L:\grussmay\20180627';   % root data folder 
cpath = cd;             % current path
cd(rpath)       
[fname,pname] = uigetfile({'*.tif;*.tiff','TIFF'});
cd(cpath)

%%
ind = strfind(fname,'.');
fext = fname(ind(end)+1:end); % file extension

% load all found sofi data
data = [];
tfname = fname;
sofiID = ind(end)-1;
sofiType = fname(ind(end)-2);
if strcmp(sofiType,'i')
    type = 'Sofi ';
    dataName = fname(1:ind(end)-7);
elseif strcmp(sofiType,'c')
    type = 'Sofi c ';
    dataName = fname(1:ind(end)-8);
else
    type = 'Sofi lin';
    dataName = fname(1:ind(end)-10);
end

for k = 1:4
    tfname(sofiID) = num2str(k);
    try 
        data(k).im = imread([pname,tfname]);
    end
end

maxOrder = length(data);
% zeropad and upscale sofi images
for k = maxOrder:-1:1
    pad = 3*k + k-1;
   temp =  padarray(data(k).im,[floor(pad/2) floor(pad/2)],0,'pre');
   temp =  padarray(temp,[ceil(pad/2) ceil(pad/2)],mean(temp(:)),'post');
   % crop the data to have only square image
   minSize = min(size(temp));
   temp = temp(3*k+1:minSize-3*k,3*k+1:minSize-3*k);
   data(k).im = temp;
   data(k).imUp = imresize(temp,size(data(end).im));
   data(k).name = [type,num2str(k)];
   data(k).fname = dataName;
end
lgd = [];
for k = 1:length(data)
    lgd{k} = ['SOFI ',num2str(k)];
end
% EXPERIMENTAL PARAMETERS
NA = 1.49; 
lam = 0.58; % wavelength [um]
dx = 0.056;  % pixel size [um]

k0 = 2*pi/lam;
kmax = pi/dx;
kc = sqrt(2)*NA*k0/kmax; % iSIM sqrt(2) resolution improvement
% kc = NA*k0/kmax;

%% plot all SOFI images

% dynamic plot of SOFI data (1 to 4th order)
% data is an array of SOFI struct
plotSOFI(data)

%% pick line in SOFI 2 and do line plots
figure(234324);imagesc(data(2).imUp)
line = imline(gca);
t = get(line,'children');
p1 = [t(1).XData ,t(1).YData];
p2 = [t(2).XData ,t(2).YData];
linewidth = 3;

figure; hold on
prof = []; lgd = [];
for k = 1:length(data)

    data(k).prof = getLineProfile(data(k).imUp,p1',p2',3);

    data(k).profX = linspace(0,1,length(data(k).prof));
    plot(data(k).profX,data(k).prof)
    
end
hold off;
legend(lgd)

%% pick line and fit a gaussian
figure(234324);imagesc(data(2).imUp)

nLines = 5; temp = [];
for k = 1:nLines
    line = imline(gca);
    temp{k} = get(line,'children');
end

%

figure; hold on
prof = []; lgd = []; res=[]; res2=[];
for n = 1:nLines
    t = temp{n};
    for k = 1:length(data)
    
        p1 = [t(1).XData ,t(1).YData];
        p2 = [t(2).XData ,t(2).YData];
        linewidth = 3;
        data(k).prof = getLineProfile(data(k).imUp,p1',p2',3);
    
        data(k).profX = linspace(0,1,length(data(k).prof));
         plot(1:length(data(k).prof),data(k).prof)
        [res(k),f,x] = fitGauss(data(k).prof);
        res(k,n) = res(k).*dx/k;
        if (res(k,n) > 0.5)
            res(k,n) = nan;
        end
        res2(k,n) = getFWHM(data(k).prof).*dx/k;
        x = 1:length(data(k).prof);
%          plot(x,f.a*exp(-((x-f.x0).^2./(2*f.s^2)))+f.b,'--','linewidth',2)
    end
end
hold off;
% legend(lgd{2:end})

figure;errorbar(1:4,mean(res2,2),std(res2,[],2),'linewidth',2); hold on;
plot(1:4,0.5*lam./(sqrt(2)*(1:4)*NA))
legend('Measured FWHM','Theoretical resolution')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMAGE PROCESSING : THIS SECTION IS DEPRECATED
% NEW ALGORITHM PRESENT IN MEASURERESOLUTION GIT

% dcorrelation measure
% requires the measureresolution git
addpath('..\measureresolution\funcs_ad');

for k = 1:length(data)
    [data(k).d,data(k).r] = estimDcor(data(k).im,10+k);
end

%% resolution estimate
g = linspace(0.001,10,30); d = [];
order = 2;
temp = data(order).im;

for k = 1:length(g)
   d(:,k) = estimDcor(temp-imgaussfilt(temp,g(k)));
end

 figure ; plot(d); legend(num2str(g'))

%% 
SNR = []; res = [];
for k = 1:size(d,2)
    [ma,ind] = max(d(:,k));
    if ind < length(d(:,k))/2
        res(k) = r(ind);
        SNR(k) = ma;
    else
        res(k) = 0;
        SNR(k) = 0;
    end
end

[q,ind] = max(SNR);
resol = res(ind)
figure;plot(res);hold on ;plot(SNR);  hold off
figure; plot(r,d);hold on ;plot([resol resol],[0 1]); hold off
%% threshold free resolution estimate based on decorrelation curves
im = data(3).im;
r0 = linspace(0,0.5,20);
x = linspace(-1,1,size(im,1));
[X,Y] = meshgrid(x);
m = zeros(length(r0),1); ind = ones(length(r0),1);
dcor = zeros(50,length(r0));

 for k = 1:length(r0)
     mask = X.^2 + Y.^2 > r0(k).^2;
    dcor(:,k) = estimDcor(real(ifftshift(ifftn(ifftshift(mask.*fftshift(fftn(fftshift(im))))))));
    [maxval,pos] = findpeaks(dcor(:,k));
    if ~isempty(maxval); [m(k),maxind] = max(maxval); end
    if ~isempty(pos); ind(k) = pos(maxind); end
 end
 
 r = linspace(0,1,50);
 figure;plot(r0,m);hold on;plot(r0,r(ind)); hold off
kc = max(r(ind));
th = linspace(0,2*pi,100);
xkc = kc.*sin(th); ykc = kc.*cos(th);
figure;imagesc(x,x,log(abs(fftshift(fftn(fftshift(im))))+1));hold on; plot(xkc,ykc,'w'); hold off
%% Filter according to kc
figure(345);imagesc(im)
mask = X.^2 + Y.^2 < kc^2;
imf = real(ifftshift(ifftn(ifftshift(mask.*fftshift(fftn(fftshift(im)))))));
figure(346);imagesc(imf)
%% compare Fourier radial profile with dcor curves
col = {'r','g','b','k'};
figure(314531); plot([0 0],[0 0]);hold on
for k = 1:4
%     temp = getRadAvg(log(abs(fftshift(fftn(fftshift(data(k).im))))+1));
%     data(k).radAvg = (temp - min(temp))./(max(temp)-min(temp));
    r = linspace(0,k,length(data(k).d));
    r2 = linspace(0,k,length(data(k).r));
%     plot(r,data(k).d)
    plot(r,data(k).d,'color',col{k},'linewidth',2)
%     plot(r,data(k).fit,'--')
     plot(r2,data(k).r,'--','color',col{k},'linewidth',2)
end
hold off

%% back fitting of 0.5 of data(k).r

for k = 1:4
    
    r = linspace(0,k,length(data(k).r));
    m = data(k).r(end/2)./(k/2);
    b = data(k).r(end/2);
    data(k).fit = 2*b - m.*r;
    
end

%%
figure(34534); hold on
for k = 1:length(data)
    r = linspace(0,k,length(data(k).d));
    [a,ind] = max(data(k).d);
    
    plot(r,data(k).d,'linewidth',2);
    cmax(k) = r(ind);
    SNR(k) = a/data(k).d(end);
    
end; hold off
% legend(lgd);
hold on
for k = 1:length(data)
    [a,ind] = max(data(k).d);
%     plot([0,max(r)],[(a+data(k).d(end))/2 (a+data(k).d(end))/2]);
    
end

hold off
figure(2532);
plot(1:length(data),dx./((1:length(data)).*res)); hold on; plot(1:length(data),log(SNR)); hold off ; legend('resolution','SNR')
%% Fourier filtering
kc = 0.5;
figure(11)
Nc = length(data);
for k = 1:length(data)
    [X,Y] = meshgrid(linspace(-k,k,size(data(k).imUp,1)));
    mask = (X.^2 + Y.^2 < (k.*kc)^2);
    data(k).imf = real(ifftshift(ifftn(ifftshift(mask.*fftshift(fftn(fftshift(data(k).imUp)))))));
    subplot(ceil(Nc/ceil(sqrt(Nc))),ceil(sqrt(Nc)),k);imagesc(data(k).imf);colorbar
end

%% Fourier analysis

figure;hold on;

radAvg = [];
for k = 1:length(data)
    data(k).fft = (fftshift(fftn(fftshift(data(k).imUp))))./numel(data(k).imUp);
    data(k).radAvg = getRadAvg(log(abs(data(k).fft)+1));
    plot(linspace(0,k,length(data(k).radAvg)),data(k).radAvg)
end

xlabel('Normalized frequencies')
NA = 1.49; lam = 0.58; k0 = sqrt(2)*2*pi/lam; 
dx = 0.056; kmax = pi/dx;
kc = NA*k0/kmax;
for k = 1:4
    plot(k.*[kc kc],[min(data(k).radAvg) max(data(1).radAvg)])
end
hold off

%% Noise analysis

figure; hold on
co = {'b','g','k','m'};
for k = 1:4
    temp = single(data(k).imUp);
    % estimate of noise & background
    noise = findNoise(temp,100);
    
    % compute spatial mask dividing noise/bckr pixels and signal pixels
    mask = temp > mean(noise(:)) + std(noise(:));
    mask = imresize(imresize(mask,0.25,'nearest'),size(temp),'bilinear');
    data(k).mask = mask;
    
    % radial circular average of fft
    S = fftshift(fftn(fftshift(temp.*mask)))./(numel(temp));
    S = log(abs(S)+1);
    N = fftshift(fftn(fftshift(temp.*(1-mask))))./(numel(temp));
    N = log(abs(N)+1);
    rS = getRadAvg(S); rN = getRadAvg(N); 
    data(k).radAvg = rS; data(k).radAvgN = rN;
    
    % crossing of curves defines achieved resolution ?
    plot(linspace(0,k,length(rS)),rS,'linewidth',2)
    plot(linspace(0,k,length(rN)),rN,'g--','linewidth',2)
end

% legend(lgd)
xlabel('Normalized frequencies')
NA = 1.49; lam = 0.58; k0 = 2*pi/lam; 
dx = 0.056; kmax = pi/dx;
kc = sqrt(2)*NA*k0/kmax;
for k = 1:4
    plot(k.*[kc kc],[min(data(k).radAvg) max(data(k).radAvg)])
end
hold off

%% Fourier filtering
figure
for k = 2
    temp = single(data(k).im(51:end-50,51:end-50));
    [X,Y] = meshgrid(linspace(-1,1,size(temp,1)));
    mask = X.^2 + Y.^2 < (0.4*k*kc)^2;
    ftemp = real(ifftshift(ifftn(ifftshift(mask.*fftshift(fftn(fftshift(temp)))))));
    subplot(121);imagesc(temp);colorbar
    subplot(122);imagesc(ftemp);colorbar
end
 %% plot all SOFI-FFT images

figure
Nc = length(data);
for k = 1:length(data)
    subplot(ceil(Nc/ceil(sqrt(Nc))),ceil(sqrt(Nc)),k);imagesc(log(abs(data(k).fft)+1));colorbar
    title(['SOFI ',num2str(k)])
end
