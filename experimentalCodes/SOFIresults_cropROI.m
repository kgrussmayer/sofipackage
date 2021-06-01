clear all;
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20140328_OnTimeSeries\1\GreenND02_405ND05_365ND09_DU897_BV_1557_wsize10_fwhm2bl0dc0_start101_sub1100';
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20140328_OnTimeSeries\1\GreenND02_405ND05_365ND09_DU897_BV_1557_wsize100_fwhm2bl0dc0_start101_sub1100';
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20140328_OnTimeSeries\1\GreenND02_405ND05_365ND09_DU897_BV_1557_wsize200_fwhm2bl0dc0_start101_sub1100';
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20140328_OnTimeSeries\1\GreenND02_405ND05_365ND09_DU897_BV_1557_wsize500_fwhm2bl0dc0_start101_sub1100';
settings.io.imagePath = 'E:\tlukes\data_SOFI\20140328_OnTimeSeries\1\GreenND02_405ND05_365ND09_DU897_BV_1557_wsize1000_fwhm2bl0dc0_start101_sub1100';

% settings.io.imagePath = 'E:\tlukes\data_SOFI\20111209_Alexa647_TIRF\001_Luc247_MONO_0548_wsize10_fwhm2bl0dc0_start101_sub1100';
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20111209_Alexa647_TIRF\001_Luc247_MONO_0548_wsize100_fwhm2bl0dc0_start101_sub1100'; 
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20111209_Alexa647_TIRF\001_Luc247_MONO_0548_wsize200_fwhm2bl0dc0_start101_sub1100';
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20111209_Alexa647_TIRF\001_Luc247_MONO_0548_wsize500_fwhm2bl0dc0_start101_sub1100';
% settings.io.imagePath = 'E:\tlukes\data_SOFI\20111209_Alexa647_TIRF\001_Luc247_MONO_0548_wsize1000_fwhm2bl0dc0_start101_sub1100';
% str = '*.mat';
str = '*.tif';

% fnames = {'001_DU897_BV_0674'}; % enter file name manually
fnames = getnamesdir(settings.io.imagePath,str)'; % search for the file names automatically
  
roix = [516,516, 433,433, 393,393,543,543,481,481];
roiy = [306,306,260,260,266,266,296,296,225,225];

roix = [184,184];
roiy = [382,382];

wsize = 200;
count = 1;

outputPath = [settings.io.imagePath,filesep,'rois'];

if (~exist(outputPath,'dir'))
    mkdir(outputPath);
end

c1=now;
for ii = 1:numel(fnames)
    
    disp(['Processing file number: ',num2str(ii),' from ',num2str(numel(fnames))]);
    if strcmp(fnames{ii}(end),'3') %crop ROI for 2nd order
        im = imread([settings.io.imagePath,filesep,fnames{ii},'.tif']);
        imcrop = im(roiy(count)+1:roiy(count)+wsize,roix(count)+1:roix(count)+wsize);
        imwrite(imcrop,[outputPath,filesep, fnames{ii},'_roi.tif']);
        count = count+1;
    end
    
end

c2=now;
% set(0, 'DefaultFigureVisible', 'on');
disp(['Total time [hours]:',num2str((c2-c1)*24)])