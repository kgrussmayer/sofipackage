% SOFI toolbox config file
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2016a
%
% config file contains path and file names of images to be processed as
% well as all the image processing settings

addpath(genpath('utils'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT FOLDERS AND FILENAMES %% CALIBRATION FILE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.io.imagePath = 'data\sofi3d\test_main_dat\20210620\input\data';
settings.io.imageName = []; % if empty, all the files in the path will be processed
settings.io.pnc = 'data\sofi3d\test_main_dat\20210620\input\calibration'; % calibration
% input file format - the same for all image inputs
% '.tif'; '.tiff' ; '.bin' ; '.dat'
settings.io.fext = '.dat';

select = []; % '*_ND04_582_75_50ms*' if select is empty, the program will read automaticly all the measurement names in the path
                % if select is a *string* enclosed by wildcarts, only folders with names containing
                % the select will be selected
                
%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%
settings.cal.logsize = 2.1; % size of the Laplacian of Gaussian filter
settings.cal.alol = 15;     % lower limit of the segments area
settings.cal.aupl = 100;    % upper limit of the segments area
settings.cal.maxshiftx = 200; % maximum coregistration shift in pixels

settings.cal.figs = 1;      % export calibration figures {0,1}
settings.cal.roix = [];     % roi used in calibration procedure
settings.cal.roiy = [];     % roi used in calibration procedure
settings.cal.px_tol = 15;	% tolerance for tentative correspondences
settings.cal.bgth = -20;    % threshold, for backgroudn suppression
settings.cal.order = 2;     % cumulant order used for registration
settings.cal.saveCal = 0;
settings.cal.correctForPSFshape = 0; % apply an extra factor to the calibration to take into account the shape of the psf


%%%%%%%%%%%%%%%%%%%%%%%%
% CUMULANT CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%
settings.sys.orders = 1:3; % sofi orders to be calculated (cumulant orders)
settings.sys.subseqlength = 1000; % length of a subsequence used for cumulant calculation
settings.sys.start = 100; % start from frame start; if empty, it starts from the first frame
settings.sys.end = []; % end at frame end; if empty, it reads until the last frame

%--- SOFI 3D ---%
settings.sys.subseqstep = 1000; % size of a sliding window, if subseqstep = subseqlength => no overlap, if subseqstep < subseqlength => overlap of subsequences
settings.sys.nplanes = 8; % number of planes in the setup
settings.sys.nss = []; % first n subsequences to be processed, if empty = process all
settings.sys.fext = settings.io.fext; % '.tif'; '.tiff' ; '.bin' ; '.dat'
% channel weights []
settings.sys.ch_weights = [0.93, 0.88, 0.80, 1.00; 0.84, 0.81, 0.85, 0.99];
% for more options see funcs/sofi/reorderData.m
% orange [166.6, 151.2, 139.2, 154.5;141.75, 147.2, 148.1, 181]
% red established Dec. 2016 average of 6 scans [0.88, 0.93, 0.73, 1; 0.88, 0.82, 0.80, 0.94]
% orange established Dec. 2016 average of 4 scans [0.93, 0.88, 0.80, 1.00; 0.84, 0.81, 0.85, 0.99]
% the channel weights for the prism setup should be [2, 4, 6, 8; 7, 5, 3, 1]
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST PROCESSING SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fwhm in [x, y, z], in pixels, [3 3 2] is default
settings.dec.fwhm = [3 3 2];

% Richardson Lucy deconvolution Matlab implementation as standard option
settings.dec.iter = 10; % number of iterations for Richardson-Lucy deconvolution, 5 is default
settings.dec.lin = 1; % turn the linearization step on/off (on = 1, off = 0)
settings.dec.medfilt = 1; % apply median filtering (on = 1, off = 0)
settings.dec.avgn = 1; % apply averaging before deconvolution (1 = no averaging, 2 = average two consecutive SR images etc.)
settings.dec.orders = settings.sys.orders; % orders to be deconvolved
settings.dec.nplanes = settings.sys.nplanes;
settings.dec.reconvolve = 0;
settings.dec.psfmodel = 'gaussian'; % 'gaussian', 'airy'

%--- 3D SOFI ---%
settings.dec.cropData = 0; % crop the data to usable FOV based on coregistration
settings.dec.axialMirror = 1; % axially mirror the stack prior to deconvolution
settings.dec.driftCorr = 1; % apply drift correction based on sofi2 data prior to deconvolution
settings.dec.estimRes = 1; % compute resolution estimate provided the ImDecorr git is installed
settings.dec.useEstimate = 0; % use resolution estimate to set deconvolution fwhm
%---------------%

%%%%%%%%%%%%%%%%%%%
% OUTPUT SETTINGS %
%%%%%%%%%%%%%%%%%%%
settings.io.bits = 16;          % number of bits of the output tif file {8,16}
settings.io.figs = 1;           % create figures yes/no {1,0}
settings.io.figshow = 0;        % show figures yes/no {1,0}
settings.io.figsave = 1;        % save figures (which figures - please explain difference in scaling) yes/no {1,0}
settings.io.figformat = 'png';  % figure will be saved in the specified format {'fig','png','pdf'}
settings.io.matsave = 0;        % save results (what results?) in matfile (can take a lot of memory on the drive)
settings.io.seq = 0; % first n 3D super-resolved images that will be saved in a sequence
settings.io.nowID = datestr(now,'yyyymmdd_HHMM');

% Output path where results will be saved
outputpath = [settings.io.imagePath];

% standard outputpath
% settings.io.outputpath = [outputpath,filesep,...
%     '_start',num2str(settings.sys.start),...
%     '_subseqlength',num2str(settings.sys.subseqlength),...
%     '_step',num2str(settings.sys.subseqstep),...
%     '_nss',num2str(settings.sys.nss),...
%     '_fwhmx',num2str(settings.dec.fwhm(1)),...
%     'y',num2str(settings.dec.fwhm(2)), ...
%     'z',num2str(settings.dec.fwhm(3)),...
%     '_iter',num2str(settings.dec.iter), ...
%     '_recon',num2str(settings.dec.reconvolve)];

% test outputpath
settings.io.outputpath = 'data\sofi3d\test_main_dat\20210620\results';

% RUN EXPERIMENT
run_sofi3d;
