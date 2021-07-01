function N = loadbinfileLength(pn,fn)
% load stack of images from binary file, function assumes that binary
% images are in the same folder as "setting.mat" which contains information
% about acquisition
% Tomas Lukes, tomas.lukes@epfl.ch

info=dir([pn fn]);

load([pn 'settings.mat']);
nx=cam1.ROIPosition(3);
ny=cam1.ROIPosition(4);

N=info.bytes/2/nx/ny;
% eof