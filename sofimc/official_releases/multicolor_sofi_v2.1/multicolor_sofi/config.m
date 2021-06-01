% Configuration file for Multicolor SOFI
%
% I/O settings 
%
% Edit path and file names to data and calibration, if the calibration data
% exists. The algorithm expects two input files (each of them supposed to be 
% a tiff file containing a sequence of images from camera 1 - first files, 
% camera 2 - second file).

% Path to image files (tiff files)
settings.io.pn = 'test_data';
% first file name (camera 1)
settings.io.fn1 = 'Alexa488_Atto565_Alexa647_ZT594RDC_Ion400-thin_PB_80s_reflected_001';
% second file name (camera 2)
settings.io.fn2 = 'Alexa488_Atto565_Alexa647_ZT594RDC_Ion400-thin_PB_80s_transmitted_001';
% input file format
settings.io.fileExtension = '.tif'; % expected '.tiff' or '.tif'

% Path to calibration files 
settings.io.pnc = 'test_data\calibration';

% Calibration file names
% first calibration file name (camera 1)
settings.io.fnc1 = 'calibration_reflected';
% second calibration file name (camera 2)
settings.io.fnc2 = 'calibration_transmitted';

% Output path where results will be saved
settings.io.outputpath = [settings.io.pn];

% Simulated ("simul") or experimental ("exper") data to be analyzed
% In the experimental setup, images from camera 2 will be flipped.
settings.io.dataType = 'simul'; % {'simul', 'exper'}

settings.io.bits = 16; % number of bits of the output tif file {8,16}
settings.io.figshow = 1;    % show figures yes/no {1,0}
settings.io.figsave = 1 ;   % save figures yes/no {1,0};

% Multicolor SOFI settings

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

% SOFI calculation settings
settings.sys.orders = 1:2; % sofi orders to be calculated (cumulant orders)
settings.sys.wsize = 500; % number of frames for the subsequence (input 
% image sequence is processed in sub-sequences / batches)
settings.sys.sub = []; % evaluate only first n frames (for quick preview)
settings.sys.start = 1; % start frame
% settings.sys.end = []; % if empty use the whole sequence
settings.sys.flattening = 'after';  % flattening {'before', 'after'} spectral unmixing 

% Calibration settings
settings.cal.order = 2; % cumulant order used for registration
if strcmp(settings.io.dataType,'exper')
    % threshold, for backgroudn suppression
    % real samples => -20 % simulations =>1 
    settings.cal.bgth=-20; 
else
    settings.cal.bgth=1;
end
settings.cal.beads = 1; % 1 use beads to coregister; 0 use data to coregister
settings.cal.logsize = 2.1; % size of the Laplacian of Gaussian filter
settings.cal.alol = 15;  % lower limit of the segments area
settings.cal.aupl = 200; % 200 upper limit of the segments area
settings.cal.maxshiftx = 200; % maximum coregistration shift in pixels

settings.cal.figs = 1; % export calibration figures {0,1}
settings.cal.roix = []; % roi used in calibration procedure
settings.cal.roiy = []; % roi used in calibration procedure
settings.cal.px_tol=15; % tolerance for tentative correspondences

% Post processing settings
settings.dec.fwhm = 2.3; % fwhm in x/y 
settings.dec.iter = 20; % number of iterations for deconvolution
settings.dec.lin = 1; % turn the linearization step on/off (on = 1, off = 0)
settings.sys.lincoeff = [1, 1/2, 1/3, 1/4];  % first, second... n-th order
settings.deconv_method = 'lucy'; % {'lucy', 'bsofi'}