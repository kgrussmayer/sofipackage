%% unmixing testing

% use this file to open already processed MC-sofi data and recompute the
% unmixing => 

% load the mc-sofi data
addpath('funcs')
% The data to be unmixed should be saved as a 3D stack
[fname,pname] = uigetfile({'*.tif;*.tiff','TIFF'});

data = [];
for k = 1:3
    data(:,:,k) = imread([pname,fname],k);
end

disp([pname,fname])

%% Spectral unmixing - using 2nd order cumulant - 3 colors
disp('Spectral unmixing')
% assign the transmission coefficient T (proportion in transmission channel)
T1 = 0.02;
T2 = 0.35;
T3 = 0.98;
% T coefficients first theo/exp
% Alexa488 0.02/0.03
% JF549 0.16/.37
% AbberiorFlip565 0.26/-
% Atto565 0.35/0.44 
% Alexa568 0.47/0.57
% Alexa647 0.98/0.99

% reimplement the option to read coefficients from file!
% % get R1,R2,R3 from calibration and filename information
% if strcmp(settings.io.dataType,'exp')
%     [R1,R2,R3] = getUnmixCoef(settings.io.pnc,settings.io.fn1);
% else
%     load([settings.io.pn,filesep,'unmixCoef_',rList(k).id]) ; % load R1 R2 and R3 into the workspace
% end

%     R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1]; % old formulation
R = [T1^2 T2^2 T3^2; 
    T1*(1-T1) T2*(1-T2) T3*(1-T3); 
    (1-T1)^2 (1-T2)^2 (1-T3)^2];
Rinv=R^(-1);

umx1c2=(Rinv(1,1)*data(:,:,3)+Rinv(1,2)*data(:,:,2)+Rinv(1,3)*data(:,:,1));
umx2c2=(Rinv(2,1)*data(:,:,3)+Rinv(2,2)*data(:,:,2)+Rinv(2,3)*data(:,:,1));
umx3c2=(Rinv(3,1)*data(:,:,3)+Rinv(3,2)*data(:,:,2)+Rinv(3,3)*data(:,:,1));


%% Flattening after unmixing if needed
    % prepare to also save the unflattened, unmixed second order data
    im_rgb_sofi2NoFlat = mergeToRgb(mean(umx3c2,3), mean(umx2c2,3), mean(umx1c2,3));
    % flattening after unmixing
    umx_c2 = {umx1c2, umx2c2, umx3c2};
    umx_c2_flat = sofiAllFlattenMC(umx_c2,2);
    umx1c2 = umx_c2_flat{1};
    umx2c2 = umx_c2_flat{2};
    umx3c2 = umx_c2_flat{3};

%% Merging unmixed cumulants into RGB image
im_rgb_sofi2 = mergeToRgb(mean(umx3c2,3), mean(umx2c2,3), mean(umx1c2,3));

figure(4);
image(im_rgb_sofi2);
title('Unmixed second order cumulants');
set(gca,'xtick',[],'ytick',[]);

r = im_rgb_sofi2(:,:,1);
g = im_rgb_sofi2(:,:,2);
b = im_rgb_sofi2(:,:,3);

figure(6)
subplot(131);
rf = zeros(size(r,1),size(r,2),3) ; rf(:,:,1) = r;  rf(:,:,3) = r; 
imshow(rf./max(rf(:)));
title('Red channel');
subplot(132);
gf = zeros(size(g,1),size(g,2),3) ; gf(:,:,2) = g;  gf(:,:,1) = g; 
imshow(gf./max(gf(:)));
title('Green channel');
subplot(133);
bf = zeros(size(b,1),size(b,2),3) ; bf(:,:,3) = b;  bf(:,:,2) = b; 
imshow(bf./max(bf(:)));
title('Blue channel');
set(gcf,'position',[465   555   1200   350])
    
