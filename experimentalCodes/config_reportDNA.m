% 2D SOFI - batch process all input files
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2015a
%
% config file contains path and file names of images to be processed as
% well as all the image processing settings

clear all;
% close all;clc;

%
%%% SPECIFY DATA PATH, SPECIFY NAME OF THE FILE
%

% path to data files - folder with all the measurements of one sample
% imagePath = 'C:\Users\Lukestom\Documents\PhD\SOFI method\SOFI Measurements\20121207_C2C12_ABTOMM20_Al647_formfixed';
imagePath = 'E:\tlukes\data_SOFI\20121207_C2C12_ABTOMM20_Al647_formfixed';
% imagePath = 'E:\tlukes\data_SOFI\20111209_Alexa647_TIRF';
% imagePath = 'E:\tlukes\data_DNA sequencing\data_DNA_Leuven_13042016';
% imagePath = 'E:\tlukes\data_DNA sequencing';

% str = '*.mat';
str = '*.tif';

% fnames = {'001_DU897_BV_0674'}; % enter file name manually
fnames = getnamesdir(imagePath,str)'; % search for the file names automatically
% fnames = {'20160127-l6ta488-00c2c3c4'};

%
%%% CHECK/ADJUST SETTINGS
%

% Preprocessing settings
settings.io.ro = 0; % reorder Nikon data (DNA data from Leuven)
settings.io.figs =1; % create figures yes/no {1,0}
settings.io.figformat ='png'; % figure will be saved in the specified format {'fig','png','pdf'}
settings.io.figshow =0; % show figures yes/no {1,0}
settings.io.figsave =1; % save figures yes/no {1,0}

% Bleaching correction settings
settings.io.blcor = 1; % bleaching correction off/on {0,1}
settings.blcor.type = 'monoexp'; % {monoexp, iir}

% Drif correction settings
settings.io.dcor = 0;% turn on or of drift correction, if on - specify path to calib file
settings.dcor.tag = '_drift'; % additional tag for drif corr. file {_drift_corr, drift}
settings.blcor.MaxCorrSamp = 500;

% Cumulant calculation settings
settings.sys.orders = 1:4;
settings.sys.wsize = 500;
settings.sys.sub = [];%[20000]; % evaluate only first n frames (for quick preview)
settings.sys.start = 101;
settings.sys.jk = 0; % turn on/off the Jacknife SNR estimation
settings.sys.block = 1;
settings.sys.pxy = 96.0384; % projected pixel size (in xy) [nm]  sofisetup = 96.0384, Hendriksetup = 104.8

% Deconvolution settings
settings.dec.fwhm = 2;
settings.dec.iter = 5;

% Molecular parameters
settings.molpar.thresh = 0.2;
settings.molpar.run = 1;

% FRC calculation
settings.frc.run = 0;
settings.frc.orders = 3;
settings.frc.bcgsub = 1.3;
settings.frc.pixelsize = settings.sys.pxy; % projected pixel size (in xy) [nm]  sofisetup = 96.0384, Hendriksetup = 104.8

% Estimate Ton
settings.ton.run = 0;
settings.ton.wsize = 100;
settings.ton.numtau = 15;

% Summary report
settings.rep.run = 1;

% I/O settings 
mname = fnames{1};
settings.io.outputpath = [imagePath,filesep,mname,'_wsize',num2str(settings.sys.wsize),'_fwhm',num2str(settings.dec.fwhm),'bl',...
    num2str(settings.io.blcor ),'dc',num2str(settings.io.dcor)];% output folder for results
if ~isempty(settings.sys.sub)
settings.io.outputpath = [settings.io.outputpath,'_sub',num2str(settings.sys.sub)];    
end
settings.io.bits = 16; % number of bits of the output tif file {8,16}
% settings.io.roi = [201:300;201:300]; % region of interest to load - keep empty [] if the whole image should be used
settings.io.roi = []; % ROI to be loaded, drif and bleach corrected
settings.io.roisx = {}; % ROI (image colums) to be processed by SOFI (example "settings.io.roisx = {61:360}")
settings.io.roisy = {}; % ROI (image rows) to be processed by SOFI (example "settings.io.roisy = {61:360}")
