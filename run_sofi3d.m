% 3D super-resolution optical fluctuation imaging (SOFI)
% author: Tomas Lukes, tomas.lukes@epfl.ch
% Tested MATLAB version 2016b 
%
% Load all the processing settings from config file  
% (see experiments)

% note: This script run a batch processing for all 3D SOFI data specified in the batch 
% 1) 3D raw cumulants of specified order are calculated
% 2) 3D raw cumulants are deconvolved, linearized
% 3) Results are save into automatically generated folder 

addpath(genpath('utils'));
% Reconstruction method
% 'sofi2d'
% 'sofi3d'
% 'sofimc' - Multicolor
% 'sofibp' - Biplane
settings.method = 'sofi3d';

fnames = settings.io.imageName; % if fnames is empty, all files in the path will be processed

%% Initiate

if isempty(fnames) && strcmp(settings.io.fext,'.bin')
    fnames = getnamesdirPattern(settings.io.imagePath,select); % search for the file names automatically
end
if ~strcmp(settings.io.fext,'.bin')
    temp = dir(settings.io.imagePath);
    temp(1:2) = [];
    temp([temp(:).isdir]) = [];
    fnames = [];
    for k = numel(temp):-1:1
        if ~isempty(strfind(temp(k).name,'.info'))
            temp(k) = [];
        end
    end
    for k = 1:numel(temp)
        fnames{k} = temp(k).name;
    end
end

% name of the output file - generate automatically TODO
folderName = strsplit(settings.io.imagePath,filesep);
if ~isempty(folderName{end})
    folderName = folderName{end};
else
    folderName = folderName{end-1}; 
end

numFrames = numel(fnames);
settings.io.folderName = folderName; 

% general data names
if strcmp(settings.io.fext,'.bin')
    settings.io.fn1 = 'data1.bin';
    settings.io.fn2 = 'data2.bin';
    settings.io.fn3 = 'data1.bin'; 
else
    settings.io.fn1 = [fnames{1}];
end

for k = 1:numel(fnames)
   disp(['File # ',num2str(k),', ',fnames{k}]) 
end

%% Run multiplane 3D SOFI

c1 = now;
for frameNumber = 1:numFrames

    settings.io.fileName = fnames{frameNumber};
    
    %%% Calculate 3D cumulants
    disp(['Processing file # ',num2str(frameNumber)])
    [sofizt,settings.sys,settings.cal] = sofi3d_process(...
        settings.sys,settings.cal,settings.io,...
        [settings.io.imagePath,filesep,fnames{frameNumber}]);
    
    c2 = now;
    %%% Deconvolution and linearization
    [sofizt,sofid,settings.dec] = sofi3d_postproc(...
        sofizt,settings.dec,settings.cal,settings.io); 
end
c3 = now;

% set(0, 'DefaultFigureVisible', 'on');
disp(['The whole processing time [hours]:',num2str((c3-c1)*24)])
disp(['Post processing part [hours]:',num2str((c3-c2)*24)])

% Save settings and processing info
settings.io.proctime = (c3-c1)*24;

save([settings.io.outputpath,filesep,'config.mat'],'settings'); % write to a mat file
descfile_sofi3D(settings); % write to a text file
