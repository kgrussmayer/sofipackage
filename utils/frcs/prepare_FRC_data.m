
% im1 = mean(sofi{1}(:,:,1:2:end),3);
% im2 = mean(sofi{1}(:,:,2:2:end),3);

im1 = mean(sofi_lin{3}(:,:,1:2:end),3);
im2 = mean(sofi_lin{3}(:,:,2:2:end),3);

% bcgsub = 0.2;
% im1 = im1 - bcgsub*median(im1(:)); 
% im2 = im2 - bcgsub*median(im2(:));
% im1(im1<0) = 0;
% im2(im2<0) = 0;

showxsubs(0,im1,im2);
% save('PALM_160127_pXJ41_CD4_CS1_PSCFP2_Jurkat_fix_ser01_125x380001_testFRC.mat','im1','im2');
% save('PALM_160224_pXJ41_CD4_wt_mEOS2_Jurkat_fix_ser02001_testFRC.mat','im1','im2');

% save('PALM_focalAdhesions_paxilin_psCFP2_mov3_roi1_lin2_wsize500_fwhm3bl0dc1_ord3.mat','im1','im2');
% save('DNA_20160127-l6ta488-00c2_ord2.mat','im1','im2');

% save('20121207_C2C12_ABTOMM20_Al647_formfixed_001_DU897_BV_0674_sofi1.mat','im1','im2');
save('20111209_Alexa647_TIRF_001_Luc247_MONO_0548_wsize1000_fwhm2.2_sofi3.mat','im1','im2');