% Configuration file for Multicolor SOFI

% I/O settings 
%
% determine what kind of data is analyzed
settings.io.dataType = 'sim'; % options are 'sim' or 'exp'

% Path to calibration files 
% settings.io.pnc = 'D:\Kristin\simulation\Alexa488_Atto565_Alexa647_ZT594RDC\calibration_ZTC594';
settings.io.pnc = 'D:\Kristin\201811-simulation\Alexa488_Atto565_Alexa647_ZT594RDC_Ion400-thin\calibration';

% Calibration file names
if strcmp(settings.io.dataType,'exp')
    settings.io.fnc1 = 'reflected';
    settings.io.fnc2 = 'transmitted';
else
    settings.io.fnc1 = 'calibration_reflected';
    settings.io.fnc2 = 'calibration_transmitted';
end

if strcmp(settings.io.dataType,'exp')
    settings.io.fileExtension = '.tiff'; % options are '.tiff' or '.tif'
else 
    settings.io.fileExtension = '.tif';
end

% Images (tiff files)
% settings.io.pn = 'D:\Kristin\simulation\Alexa488_Atto565_Alexa647_ZT594RDC';
settings.io.pn = 'D:\Kristin\201811-simulation\Alexa488_Atto565_Alexa647_ZT594RDC_Ion400-thin';

% settings.io.fn1 = 'patches_Alexa488_Atto565_Alexa647_ZT594RDC_Ion400_sameParameters_reflected_000';
% settings.io.fn2 = 'patches_Alexa488_Atto565_Alexa647_ZT594RDC_Ion400_sameParameters_transmitted_000';
settings.io.fn1 = 'Alexa488_Atto565_Alexa647_ZT594RDC_Ion400-thin_PB_80s_reflected_001';
settings.io.fn2 = 'Alexa488_Atto565_Alexa647_ZT594RDC_Ion400-thin_PB_80s_transmitted_001';

% settings.io.fn1 = 'HeLa_microtubules-Alexa488_WGA-Atto565_LaminB1-Alexa647_PCA_PCD_full488nm_full561nm_0_33A_635nm_highIllInt_reflected_013';
% settings.io.fn2 = 'HeLa_microtubules-Alexa488_WGA-Atto565_LaminB1-Alexa647_PCA_PCD_full488nm_full561nm_0_33A_635nm_highIllInt_transmitted_013';

% settings.io.fn1 = 'HeLa_microtubules-Alexa488_WGA-Atto565_LaminB1-Alexa647_PCA_PCD_full488nm_full561nm_0_33A_635nm_reflected_005';
% settings.io.fn2 = 'HeLa_microtubules-Alexa488_WGA-Atto565_LaminB1-Alexa647_PCA_PCD_full488nm_full561nm_0_33A_635nm_transmitted_005';

settings.io.bits = 16; % number of bits of the output tif file {8,16}
settings.io.outputpath = [settings.io.pn];
settings.io.figshow = 0;    % show figures yes/no {1,0}
settings.io.figsave = 1 ;   % save figures yes/no {1,0}
settings.io.savestack = 0;  % save stacks yes/no {1,0}

% Multicolor SOFI settings

% only for batch processing with triplet Alexa488, X, Alexa647
% automatically determine X from filename and theo or exp coefficient version
% from exp or sim setting
settings.mc.getUnmixCoefAlexa488XAlexa647 = 0; % use automatic determination 1
% % new formulation of the matrix
settings.mc.T1 = 0.02;        
settings.mc.T2 = 0.35;        
settings.mc.T3 = 0.98;
% old formulation of the matrix ATTENTION R is not reflection coefficient
% but ratio!
settings.mc.R1 = 0.02;        
settings.mc.R2 = 0.55;        
settings.mc.R3 = 57.45;
% T coefficients first theo/exp
% Alexa488 0.02/0.03
% JF549 0.16/.37
% AbberiorFlip565 0.26/-
% Atto565 0.35/0.44 
% Alexa568 0.47/0.57
% Alexa647 0.98/0.99

% Cumulant calculation settings
settings.sys.orders = 1:3; % sofi orders to be calculated
settings.sys.wsize = 500; % number of frames for the subsequence
settings.sys.sub = [];%[20000]; % evaluate only first n frames (for quick preview)
settings.sys.start = 1; %1001
settings.sys.end = []; % if empty use the whole sequence
settings.sys.flattening = 'after';  % apply flattening before or after the spectral unmixing {'before', 'after'}

% Calibration settings
settings.cal.order = 2; % cumulant order used for registration
if strcmp(settings.io.dataType,'exp')
    settings.cal.bgth=-20; % real samples => -20 % simulations =>1 background threshold, make sure not too much noise is picked up!
else
    settings.cal.bgth=1;
end
settings.cal.beads = 1; % 1 use beads to coregister; 0 use data to coregister
settings.cal.logsize=2.1; % size of the Laplacian of Gaussian filter
settings.cal.alol=15;  % lower limit of the segments area
settings.cal.aupl=200; % 200 upper limit of the segments area
settings.cal.maxshiftx = 200; % maximum coregistration shift in pixels

settings.cal.figs = 1; % export calibration figures {0,1}
settings.cal.roix = []; % ! 1:100 if not empty it determines the roi used in calibration procedure!
settings.cal.roiy = []; % ! 1:100 if not empty it determines the roi used in calibration procedure!
settings.cal.px_tol=15; % tolerance for tentative correspondences

% Post processing settings
settings.dec.fwhm = 2.3; % fwhm in x/y 
settings.dec.iter = 20; % number of iterations for Lucy 
settings.dec.lin = 1; % turn the linearization step on/off (on = 1, off = 0)
settings.sys.lincoeff = [1, 1/2, 1/3, 1/4];  % first, second... n-th order
% note: [1/2, 1/3, 1/4] are the standard coeffiecient from theory, 
% approx. [0.75, 0.55, 0.4] are experimentally proven to be the most often
% estimated values if the adaptive linearisation is applied
settings.dec.medfilt = 1; % apply median filtering (on = 1, off = 0) 
settings.dec.avgn = 1; % apply averaging before deconvolution (1 = no averaging, 2 = average two consecutive SR images etc.)
settings.dec.orders = 2:4; % orders to be deconvolved
settings.dec.nplanes = 2;
settings.deconv_method = 'bsofi'; % {'breg_cuda', 'augLag_mb', 'lucy', 'bsofi'}