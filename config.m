% We are currently in the process of cleaning up this repository 
% and adding a proper user manual and testing - stay tuned for updates

% 2D SOFI - batch process all input files
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2016a
%
% config file contains path and file names of images to be processed as
% well as all the image processing settings

clear;
% close all;clc;

% INPUT SETTINGS
% path to data files - folder with all the measurements to process
str = '.raw';
% (.tif, .tiff, '.dat', '.raw') specify the file type by the file extension, 
% currently reads original tif(f) files up to 4Gb

settings.io.imagePath = 'insertYourPath';
settings.io.fext = str;
settings.io.W = 400;
settings.io.H = 400;
fnames = getnamesdir(settings.io.imagePath, str)'; % search for the file names automatically

settings.io.roi = []; % ROI to be loaded, drift and bleach corrected- keep empty [] if the whole image should be used
settings.io.roisx = {}; % ROI (image colums) to be processed by SOFI (example "settings.io.roisx = {61:360}")
settings.io.roisy = {}; % ROI (image rows) to be processed by SOFI (example "settings.io.roisy = {61:360}")

%%% CHECK/ADJUST IMAGE PROCESSING SETTINGS

% Preprocessing settings
settings.io.ro = 0; % reorder Nikon data (DNA data from Leuven)
settings.io.figs =1; % create figures yes/no {1,0}
settings.io.figformat = 'png'; % figure will be saved in the specified format {'fig','png','pdf'}
settings.io.figshow =0; % show figures yes/no {1,0}
settings.io.figsave =1; % save figures (which figures - please explain difference in scaling) yes/no {1,0}
settings.io.matsave =0; % save results (what results?) in matfile (can take a lot of memory on the drive)
settings.io.concatOn =0; % concatenate consecutive images yes/no {1,0} TODO: needs better explanation!
settings.io.concatSave =0; % save concatenate images in one file yes/no {1,0} TODO: needs better explanation!

% Bleaching correction settings
settings.io.blcor = 1; % bleaching correction off/on {0,1}
settings.blcor.type = 'monoexp';%'monoexp'; % {monoexp, iir}
settings.blcor.MaxCorrSamp = 5000; % TODO: needs better explanation!

% Drift correction settings TODO: please check explanation and complete it!
settings.io.dcor = 1;% turn on or off drift correction, if on 
% - specify path to drift correction file if drift correction method is other than SOFI
settings.dcor.type = 'SOFI'; % drift correction method {'TS','LBEN_PALM','SOFI'} 
% We assume the drift correction file to be : [settings.io.imageFile,settings.dcor.tag]
settings.dcor.tag = '_drift';%'driftcor'; % additional tag for drif corr. file {_drift_corr, drift,driftcor}
% TS: drift correction based on ThunderSTORM .json file (either from fiducial markers or cross-correlation)
% LBEN_PALM:
% SOFI: first cross-correlation between 100 frames subsequences (why not variable?) and first 100 frames, performed for mean of stack?? 
% then followed by cross-correlation between first and other SOFI2
% subsequences (is this the same as in biplane?)

% Cumulant calculation settings
settings.sys.orders = 1:4;
settings.sys.wsize = 1000; % substack for analysis to avoid correlations from bleaching etc
settings.sys.sub = [];%[550];%10050%[20000]; % evaluate only first n frames (for quick preview)
settings.sys.start = []; %51; start from frame start; if empty, it starts from the first frame
settings.sys.end = []; % end at frame end; if empty, it reads until the last frame 

% POST PROCESSING SETTINGS

% TODO: which part is this used for?
settings.sys.pxy = 108; % projected pixel size (in xy) [nm]  sofisetup = 96.0384, Hendriksetup = 104.8, AD-gut setup = 108

% Deconvolution/Linearization settings
settings.dec.deconv_method = 'lucy'; % {'breg_cuda', 'augLag_mb', 'lucy'}
settings.dec.fwhm = 3;
settings.dec.lin = 1; % option to turn on/off linearization
settings.dec.reconvolve = 0; % option to turn on/off reconvolution using 
% conv and the theoretical psf model "airy", or "gaussian" chosen below
settings.dec.denoise = 1; % only applies to SOFI lin

% parameters for Lucy-Richardson deconvolution
settings.dec.iter = 5;
settings.dec.psfmodel = 'gaussian';% choice of theoretical "airy", or "gaussian"

% parameters for cuda version of bregman iterative method based
% deconvolution with Gaussian noise model and Gaussian psf model
settings.bregman.iter = 10;
settings.bregman.apodize =0;
settings.bregman.lambda = 0.03;
settings.bregman.omega = 0.95;
settings.bregman.NumImages = 1;

% parameters for matlab version of augmented lagrangian based deconvolution
% with Gaussian noise model and Gaussian psf model
settings.augLag.gamma = 2000;
settings.augLag.beta = settings.augLag.gamma;
settings.augLag.alpha = 1;
settings.augLag.reltol = 1e-4;
settings.augLag.maxIter = 3;
settings.augLag.Lp = 1;
% settings.augLag.Lp = 1; % defines which Lp norm to use (hardcoded to 1 in
% the code)

% Molecular parameters
settings.molpar.run = 1; % turn on or off
settings.molpar.thresh = 0.1; % TODO: needs better explanation!
% please add an option for saving in .mat file and .tif
% use tirf yes/no as option for brightness calculation

% Estimate Ton
settings.ton.run = 1; % turn on or off
settings.ton.wsize = 300;
settings.ton.numtau = 20; % TODO: needs better explanation!

% FRC calculation
settings.frc.run = 0; % turn on or off
settings.frc.orders = 1:3;
settings.frc.bcgsub = 1;%1.3;
settings.frc.pixelsize = settings.sys.pxy; % projected pixel size (in xy) [nm]  sofisetup = 96.0384, Hendriksetup = 104.8, AD-gut setup = 108

% Jacknife SNR estimation
settings.sys.jk = 1; % turn on/off the Jacknife SNR estimation - where is it in the code?
settings.jk.orders = 1:4;
settings.sys.block = 1; % TODO: needs better explanation!

% Summary report only if bleaching correction is used
settings.rep.run = 0;

% OUTPUT SETTINGS 
mname = fnames{1};
settings.io.outputpath = settings.io.imagePath; 
%     num2str(settings.io.blcor ),'dc',num2str(settings.io.dcor),'_start',num2str(settings.sys.start)];% output folder for results
% if ~isempty(settings.sys.sub)
settings.io.outputpath = [settings.io.outputpath,'_', settings.dec.deconv_method,'_FWHM', num2str(settings.dec.fwhm), '_iter',num2str(settings.dec.iter),'_recon',num2str(settings.dec.reconvolve), '_wsize', num2str(settings.sys.wsize), ...
    '_start', num2str(settings.sys.start), '_end', num2str(settings.sys.end)];    
settings.io.bits = 16; % number of bits of the output tif file {8,16}

%TODO please describe the results/ output better here or in a readme file
% SOFI PROCESSING
% raw cumulants sofi_c
% deconvolved/linearized cumulants
% adaptive linearization sofi_lin
% standard linearization using fixed coefficients sofi

% MOLECULAR PARAMETER ESTIMATES/bSOFI
% indicate which algorithm is used - 2nd,3rd,4th order combination from bSOFI paper 
% or 2nd, 3rd order for on-time ratio as in PALM & SOFI paper
% it should be clear which flattening/linearization is used as basis for
% bSOFI at the moment it is adaptive linearization
% figures in fig folder: 
% results overview pdf: 
% .mat file:
% please clarify the units for the molecular parameter estimates - e.g. is
% it molecules/px2?
