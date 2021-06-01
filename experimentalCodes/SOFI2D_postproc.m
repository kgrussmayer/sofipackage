% 2D SOFI processing: perform flattening, deconvolution and linearization on
% precomputed raw cumulant images
% Tomas Lukes, tomas.lukes@epfl.ch 
clear all;
addpath('E:\tlukes\codesgit\myutils');
% addpath('E:\tlukes\codesold\test_cudaDeconv');
% addpath('E:\tlukes\codesold\test_cudaDeconv\linearizations');
addpath('E:\tlukes\codesold\test_deconvSOFI');
addpath('E:\tlukes\codesold\test_deconv3D_v02\funcs');
addpath('E:\tlukes\codesgit\wdt\deconvs');

config_postproc;

settings.sys.orders = 2;
settings.dec.fwhm = 2.5;
settings.dec.iter = 10;

settings.io.outputpath = [imagePath,filesep,'postproc',filesep,fnames{1},'_wsize',num2str(settings.sys.wsize),'_fwhm',num2str(settings.dec.fwhm),'bl',...
    num2str(settings.io.blcor ),'dc',num2str(settings.io.dcor),'_start',num2str(settings.sys.start),'_1'];% output folder for results

settings.io.imageName = fnames{6}; 
settings.io.imageFile = [imagePath, filesep,settings.io.imageName];
    
%
%%% Load image stack
%
% [stack,frames] = loadStack(settings);
stack = load_tifFile([settings.io.imageFile,'.tif'],settings.sys.sub,settings.io.roi);
frames=size(stack,3);
settings.sys.frames = frames;
numim = size(stack,3);
maxf = 1; %first 1000 frames
c{settings.sys.orders}=double(stack(:,:,1:maxf)); 

%
%%% Flattening
%

c=sofiAllFlatten(c,settings.sys.orders);

for ii = 1:maxf
    %         [c,fwhm]=sofiFlatten([],c,grid,settings.sys.orders);
    %         fwhms(count)=fwhm;

    for io=settings.sys.orders
        c{io}(:,:,ii)=medfilt2(c{io}(:,:,ii),[2 2]);
    %     sofi_c{io}(:,:,ii)=c{io};
    end

end


%
%%% Deconvolution
%


sofi=sofiLinearize(c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter);
% sofi_lin = sofiAdalin(c, stack, settings); 
%%
fwhm = settings.dec.fwhm;
iter = 100;
lincoeffs = [1 0.4 0.45 0.35];
fcpx = 0.95;

% [sofi]=sofiLinearize_mbm(c,fwhm,settings.sys.orders,iter,lincoeffs,fcpx,1);
% [sofi3]=sofiLinearize_mbm2(c,fwhm,settings.sys.orders,iter,lincoeffs,fcpx);
% [sofi]=sofiLinearize_mb(c,fwhm,settings.sys.orders,iter,lincoeffs,fcpx,1);
%%
% Settings for the deconvolution algorithm
settings.reltol = 1e-4;
settings.maxIter = 10;
settings.numIm = 1;
settings.lambda2 = 0.05;
settings.coeff_A = 0.1;
settings.dz = 1;

%%% mseBreg_v03

% PSF parameters
order = settings.sys.orders;

fwhm1=(fwhm)/sqrt(8*log(2));
PSF=gaussian(fwhm1*sqrt(order));
PSF=PSF(:)*PSF;
    
im = c{order};
bluredge = fspecial('gaussian',15,5); 
im = edgetaper(im,bluredge);
[imc,concur3] = mseBreg_v03(im,PSF,settings);


%%
settings.coeff_H = 1.2; % for auglag 19
settings.gamma = 0.0001; % for auglag 19
[imd2,concur4] = deconvAugLag3D_v19(im,PSF,settings);

% [im] = mseBreg_v02b(stack,PSF,settings);
% [im] = mseBreg_v02_mirroring(c{order},PSF,settings);

sofi2{order} = imc;

%%
settings.coeff_H = 1;
% settings.beta = 1;
settings.gamma = 0.01; 
settings.maxIter = 50;

[imd2,concur7] = poissBreg3(im,PSF,settings);

%%
order = settings.sys.orders;
desc.title = 'deconv comparison';
% showxsubs(desc,c{order},sofi{order},sofi2{order}.^lincoeffs(order));
temp = sofi2{order};
% temp = min(0,temp - 1*median(temp(:)));
showxsubs(desc,c{order}.^0.7,sofi{order}.^0.7,temp.^0.7,imd2.^0.7);

% showxsubs(desc,c{order},sofi2{order}.^lincoeffs(order));
sofi{order} = imd2.^0.7;
sofiSaveResults(settings,sofi);
