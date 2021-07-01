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

settings.io.imagePath = 'data\sofimc\test_main\20210615\input'; 
% first file name (camera 1)
settings.io.imageName = 'Alexa488_Atto565_Alexa647_reflected';
% second file name (camera 2) - for multicolor
settings.io.imageName2 = 'Alexa488_Atto565_Alexa647_transmitted';

settings.io.pnc = 'data\sofimc\test_main\20210615\input\calibration';
% first calibration file name (camera 1)
settings.io.fnc1 = 'calibration_reflected';
% second calibration file name (camera 2)
settings.io.fnc2 = 'calibration_transmitted';
settings.io.fext = '.tif';
% input file format - the same for all image inputs
% '.tif'; '.tiff' ; '.bin' ; '.dat'

% Simulated ("simul") or experimental ("exper") data to be analyzed
% In the experimental setup, images from camera 2 will be flipped.
% {'simul', 'exper'}
settings.io.dataType = 'simul'; % {'simul', 'exper'}

%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%
settings.cal.beads = 1;
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

% threshold, for background suppression
switch settings.io.dataType
    case 'exper' % real samples => -20
        settings.cal.bgth = -20;
    case 'simul' % simulations => 1
        settings.cal.bgth = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% CUMULANT CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%
settings.sys.orders = 1:2; % sofi orders to be calculated (cumulant orders)
settings.sys.subseqlength = 500; % length of a subsequence used for cumulant calculation
% number of frames in subsequence (input image sequence is processed in sub-sequences / batches)
settings.sys.sub = []; % evaluate only first n frames (for quick preview)
settings.sys.start = []; % start from frame start; if empty, it starts from the first frame
settings.sys.end = []; % end at frame end; if empty, it reads until the last frame
settings.sys.nplanes = 8; % number of planes in the setup
settings.sys.flattening = 'after';  % flattening {'before', 'after'} spectral unmixing 

%--- MULTICOLOR ---%
settings.sys.flattening = 'after';  % flattening {'before', 'after'} spectral unmixing
% T coefficients (theoretical / experimental)
% Alexa488 0.02/0.03
% JF549 0.16/.37
% AbberiorFlip565 0.26/-
% Atto565 0.35/0.44
% Alexa568 0.47/0.57
% Alexa647 0.98/0.99
settings.mc.T1 = 0.02;
settings.mc.T2 = 0.35;
settings.mc.T3 = 0.98;
%------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST PROCESSING SETTINGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.dec.fwhm = [2.3 2.3 3];
% Richardson Lucy deconvolution Matlab implementation as standard option
settings.dec.iter = 20; % number of iterations for Richardson-Lucy deconvolution, 5 is default
settings.dec.lin = 1; % turn the linearization step on/off (on = 1, off = 0)
settings.dec.medfilt = 1; % apply median filtering (on = 1, off = 0)
settings.dec.avgn = 5; % apply averaging before deconvolution (1 = no averaging, 2 = average two consecutive SR images etc.)
settings.dec.orders = settings.sys.orders; % orders to be deconvolved
settings.dec.nplanes = settings.sys.nplanes;
settings.dec.reconvolve = 0;
settings.dec.psfmodel = 'gaussian'; % 'gaussian', 'airy'

%--- MULTICOLOR ---%
settings.sys.lincoeff = [1, 1/2, 1/3, 1/4];  % first, second... n-th order
settings.deconv_method = 'lucy'; % {'lucy', 'bsofi'}
%------------------%

%%%%%%%%%%%%%%%%%%%
% OUTPUT SETTINGS %
%%%%%%%%%%%%%%%%%%%
settings.io.bits = 16;          % number of bits of the output tif file {8,16}
settings.io.figs = 1;           % create figures yes/no {1,0}
settings.io.figshow = 0;        % show figures yes/no {1,0}
settings.io.figsave = 1;        % save figures (which figures - please explain difference in scaling) yes/no {1,0}
settings.io.figformat = 'png';  % figure will be saved in the specified format {'fig','png','pdf'}
settings.io.matsave = 0;        % save results (what results?) in matfile (can take a lot of memory on the drive)
settings.io.seq = 0;            % first n 3D super-resolved images that will be saved in a sequence

% Output path where results will be saved

settings.io.outputpath = 'data\sofimc\test_main\20210615\result\output';

% RUN EXPERIMENT
run_sofi_multicolor;
