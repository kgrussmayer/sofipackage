% SOFI toolbox config file
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2016a
%
% config file contains path and file names of images to be processed as
% well as all the image processing settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT FOLDERS AND FILENAMES %% CALIBRATION FILE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.io.imagePath =  'D:\Data\2022\temporal\sofipackage_nobin\';
settings.io.imageName = [];
settings.io.fext = '.tif';

%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%
% input file format - the same for all image inputs
% '.tif'; '.tiff' ; '.bin' ; '.dat'

settings.io.W = 200;
settings.io.H = 200;
settings.io.roi = [];   % ROI to be loaded, drift and bleach corrected
% - keep empty [] if the whole image should be used
settings.io.roisx = {}; % ROI (image colums) to be processed by SOFI
% - example: settings.io.roisx = {61:360}
settings.io.roisy = {}; % ROI (image rows) to be processed by SOFI
% - example: settings.io.roisy = {61:360}

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/O PROCESSING SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- 2D SOFI ---%
% Bleaching correction settings
settings.io.blcor = 1;              % bleaching correction off/on {0,1}
settings.blcor.type = 'monoexp';    % {monoexp, iir}
settings.blcor.MaxCorrSamp = 5000;  % TODO: needs better explanation!

settings.io.dcor = 0;% turn on or off drift correction, if on
% - specify path to drift correction file if drift correction method is other than SOFI
settings.dcor.type = 'SOFI'; % drift correction method {'TS','LBEN_PALM','SOFI'}
% We assume the drift correction file to be : [settings.io.imageFile,settings.dcor.tag]
settings.dcor.tag = '_drift';%'driftcor'; % additional tag for drif corr. file {_drift_corr, drift,driftcor}
% SOFI: first cross-correlation between 100 frames subsequences (why not variable?) and first 100 frames, performed for mean of stack??
% then followed by cross-correlation between first and other SOFI2 subsequences
settings.io.ro = 0; % reorder Nikon data (DNA data from Leuven)
%---------------%

settings.io.concatOn = 0;   % concatenate consecutive images no/yes {0,1}
settings.io.concatSave = 0; % save concatenate images in one file no/yes {0,1}

%%%%%%%%%%%%%%%%%%%%%%%%
% CUMULANT CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%
settings.sys.orders = 1:2; % sofi orders to be calculated (cumulant orders)
settings.sys.subseqlength = 1000; % length of a subsequence used for cumulant calculation
settings.sys.start = []; % start from frame start; if empty, it starts from the first frame
settings.sys.end = []; % end at frame end; if empty, it reads until the last frame
settings.sys.sub = []; % evaluate only first n frames (for quick preview)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST PROCESSING SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fwhm in [x, y, z] in pixels, [3 3 2] is default
% multiplane setup pixel size 111 nm, interplane distance 350 nm
settings.dec.fwhm = 5;
% Richardson Lucy deconvolution Matlab implementation as standard option
settings.dec.iter = 10; % number of iterations for Richardson-Lucy deconvolution, 5 is default
settings.dec.orders = settings.sys.orders; % orders to be deconvolved
settings.dec.reconvolve = 1;
settings.dec.psfmodel = 'gaussian'; % 'gaussian', 'airy'

%--- 2D SOFI ---%
settings.sys.pxy = 156; % projected pixel size (in xy) [nm]  sofisetup = 96.0384, PALM setup = 104.8, AD-gut setup = 108
% Linearization settings
%settings.dec.lin = 1; % turn the linearization step on/off (on = 1, off = 0)
settings.dec.denoise = 1; % option to turn off/on {0, 1}

% Augmented Lagrangien deconvolution settings
settings.augLag.deconv = 0; % extra deconvolution using utils\decon\deconvAugLag3D.m
settings.augLag.gamma = 4000;
settings.augLag.FWHM = 3;
settings.augLag.beta = settings.augLag.gamma;
settings.augLag.alpha = 1;
settings.augLag.reltol = 1e-4;
settings.augLag.maxIter = 20;

% Molecular parameters
settings.molpar.run = 0; % turn on or off
settings.molpar.thresh = 0.05; % bacground suppression threshold
% Estimate Ton
settings.ton.run = 0; % turn on or off
settings.ton.subseqlength = 500;
settings.ton.numtau = 20; % TODO: needs better explanation!
% FRC calculation
settings.frc.run = 0; % turn on or off
settings.frc.orders = settings.sys.orders;
settings.frc.bcgsub = 1;
settings.frc.pixelsize = settings.sys.pxy; % projected pixel size (in xy) [nm]  
% Jacknife SNR estimation
settings.sys.jk = 0; % turn on/off the Jacknife SNR estimation
settings.jk.orders = 1:2;
settings.sys.block = 1; % TODO: needs better explanation!
% Summary report only if bleaching correction is used
settings.rep.run = 0;
%---------------%

%%%%%%%%%%%%%%%%%%%
% OUTPUT SETTINGS %
%%%%%%%%%%%%%%%%%%%
settings.io.bits = 16; % number of bits of the output tif file {8,16}

settings.io.figs = 1;           % create figures yes/no {1,0}
settings.io.figshow = 1;        % show figures yes/no {1,0}
settings.io.figsave = 1;        % save figures (which figures - please explain difference in scaling) yes/no {1,0}
settings.io.figformat = 'png';  % figure will be saved in the specified format {'fig','png','pdf'}
settings.io.matsave = 0;        % save results (what results?) in matfile (can take a lot of memory on the drive)
settings.io.seq = 0; % first n 3D super-resolved images that will be saved in a sequence

% Output path where results will be saved
settings.io.nowID = datestr(now,'yyyymmdd_HHMM');

% standard outputpath
settings.io.outputpath = [settings.io.imagePath,'output_',settings.io.nowID];
    %'_iter',num2str(settings.dec.iter),...
    %'_recon',num2str(settings.dec.reconvolve),...
    %'_subseqlength',num2str(settings.sys.subseqlength),...
    %'_start',num2str(settings.sys.start),...
    %'_end',num2str(settings.sys.end),...
    %'_',settings.io.nowID];

% RUN EXPERIMENT
run_sofi2d;
   