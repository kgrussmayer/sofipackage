% SOFI toolbox config file
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2016a
%
% config file contains path and file names of images to be processed as
% well as all the image processing settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT FOLDERS AND FILENAMES %% CALIBRATION FILE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings.io.imagePath = 'data\sofibiplane\test_main\2021618\input\data';
% if empty, image pairs in the input folder will be found automatically
settings.io.imageName = 'COS-7_Tubulin_SOFI_Flip565_biplane_reflected';
settings.io.imagePath2 = settings.io.imagePath;
settings.io.imageName2 = 'COS-7_Tubulin_SOFI_Flip565_biplane_transmitted';
settings.io.pnc = 'data\sofibiplane\test_main\2021618\input\calibration';
settings.io.fnc1 = 'Biplane_beads_calibration';
% input file format - the same for all image inputs
% '.tif'; '.tiff' ; '.bin' ; '.dat'
settings.io.fext = '.tiff';

%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%
settings.cal.beads = 0;
% 1 - use beads to coregister generates an affine transformation
%     using a stack of beads supplied as calibration and then uses SOFI 2nd
%     order cross-correlation for refinement
% 0 - uses only data to coregister, cross-correlation is calculated for the
%     SOFI 2nd order to estimate translation in x and y
settings.cal.logsize = 2.1; % size of the Laplacian of Gaussian filter
settings.cal.alol = 15;     % lower limit of the segments area
settings.cal.aupl = 100;    % upper limit of the segments area
settings.cal.maxshiftx = 200; % maximum coregistration shift in pixels

settings.cal.figs = 1;      % export calibration figures {0,1}
settings.cal.roix = [];     % roi used in calibration procedure
settings.cal.roiy = [];     % roi used in calibration procedure
settings.cal.px_tol = 15;	% tolerance for tentative correspondences
settings.cal.bgth = -20;    % threshold, for background suppression

%%%%%%%%%%%%%%%%%%%%%%%%
% CUMULANT CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%
settings.sys.orders = 1:4; % sofi orders to be calculated (cumulant orders)
settings.sys.subseqlength = 500; % length of a subsequence used for cumulant calculation
settings.sys.start = []; % start from frame start; if empty, it starts from the first frame
settings.sys.end = []; % end at frame end; if empty, it reads until the last frame
settings.sys.sub = [];

%--- BIPLANE ---%
% channel weights [reflected, transmitted]
settings.sys.ch_weights = [1, 1];
settings.sys.nplanes = 2; % number of planes in the setup
%---------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST PROCESSING SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fwhm in [x, y, z], in pixels, [3 3 2] is default
settings.dec.fwhm = [2.5 2.5 0.2];
% Richardson Lucy deconvolution Matlab implementation as standard option
settings.dec.iter = 10; % number of iterations for Richardson-Lucy deconvolution, 5 is default
settings.dec.lin = 1; % turn the linearization step on/off (on = 1, off = 0)
settings.dec.medfilt = 1; % apply median filtering (on = 1, off = 0)
settings.dec.avgn = 5; % apply averaging before deconvolution (1 = no averaging, 2 = average two consecutive SR images etc.)
settings.dec.orders = 2:4; % orders to be deconvolved
settings.dec.nplanes = settings.sys.nplanes;
settings.dec.reconvolve = 0;
settings.dec.psfmodel = 'gaussian'; % 'gaussian', 'airy'

%--- BIPLANE---%
settings.dec.cropData = 0; % crop the data to usable FOV based on coregistration
settings.dec.axialMirror = 0; % axially mirror the stack prior to deconvolution
settings.dec.driftCorr = 1; % apply drift correction based on sofi2 data prior to deconvolution
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

% Output path where results will be saved
settings.io.nowID = datestr(now,'yyyymmdd_HHMM');

% standard outputpath
% settings.io.outputpath = [settings.io.imagePath,...
%     '_subsequence',num2str(settings.sys.subseqlength), ...
%     '_iter',num2str(settings.dec.iter),...
%     '_recon',num2str(settings.dec.reconvolve),...
%     '_beads',num2str(settings.cal.beads),...
%     '_',settings.io.nowID];

% test outputpath
settings.io.outputpath = 'data\sofibiplane\test_main\2021618\result\output';

% RUN EXPERIMENT
run_sofi_biplane;
