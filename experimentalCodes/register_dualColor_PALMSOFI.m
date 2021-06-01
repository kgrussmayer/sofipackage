% load two tif files (two color channels) and PALM localizations
% estimate transformation matrix, register colorchannels

ph_mov = 'E:\tlukes\data_SOFI\20160709_PALMSOFI_dualcolor\20160511_mov2_meos2_wsize500_fwhm2bl1dc1';
ph_fox = 'E:\tlukes\data_SOFI\20160709_PALMSOFI_dualcolor\20160511_mov2_pscfp2_wsize500_fwhm2bl1dc1_start86';
ph_fox = 'E:\tlukes\data_SOFI\20160709_PALMSOFI_dualcolor\20160511_mov2_pscfp2_wsize500_fwhm2bl1dc1_start85';

fname = '20160511_mov2_meos2sofi_lin4.tif';

im_mov = double(imread([ph_mov,filesep,'20160511_mov2_meos2sofi_lin4.tif']));
im_fix = double(imread([ph_mov,filesep,'20160511_mov2_pscfp2sofi_lin4.tif']));

load('E:\tlukes\data_SOFI\20160709_PALMSOFI_dualcolor\grid_positions_meos2.mat');
load('E:\tlukes\data_SOFI\20160709_PALMSOFI_dualcolor\grid_positions_pscfp2.mat');

order = 4; % SOFI order
pxsize = 104.8/4;
figs =1;

ccf = grid_positions_meos2./pxsize;  
ccm = grid_positions_pscfp2./pxsize;

tform = fitgeotrans(ccf(:,[1 2]),ccm(:,[1 2]),'polynomial', 2);

    im_reg = imwarp(im_mov, tform,'OutputView',imref2d(size(im_fix)));
    
%     im_mov = im_mov./max(im_mov(:))*2;
%     im_fix = im_fix./max(im_fix(:))*1.5;
%     im_reg = im_reg./max(im_reg(:))*2;

if figs ==1;
    figure('Visible','On'),     
    imshowpair(im_fix, im_reg,'Scaling','joint'); %green: first image, magenta: second image
    title('registered meos2 and pscfp2 SOFI images ')  
    
    figure('Visible','On'),
    imshowpair(double(im_fix), double(im_mov),'Scaling','joint'); %green: first image, magenta: second image
    title('SOFI and PALM pts before registration')
    hold on;
    scatter(ccm(:,1),ccm(:,2),'x','MarkerEdgeColor',[1 0 1]);
    scatter(ccf(:,1),ccf(:,2),'x','MarkerEdgeColor',[0 1 0]);
    
    legend('ch1 meos2','ch 2 pscfp2');
     % save the current figure in tif file
%     fh = gcf;saveFigure(fh,outputPath,'detectedBeads','tif'); 
end

saveImageStack(im_reg ,ph_mov,[fname,'_registered'],[],16)
