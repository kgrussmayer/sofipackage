% 2D SOFI - batch process all input files
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2015a

% Load file names and processing settings from the config file
addpath(genpath('utils'));
config;
% config_local; 
 
c1=now;
warning('off','all')

if settings.io.figshow ==1
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');  
end
    
%% SOFI 2D batch 
[fnames2,~] = getFileNums(fnames); % detect numbers from input image names 
unames = uniNames(fnames,fnames2); % group unique names
settings.io.allNames = char(fnames);

for ii = 1:numel(unames)
    c(ii)=now;
    disp(['Processing file number: ',num2str(ii),' from ',num2str(numel(unames))]);
    
    if settings.io.concatOn == 1
        settings.io.imageName = unames(ii).fnames;  
    else
        if strcmp(str,'.tif') == 1 
            settings.io.imageName = unames(ii).fnames; % 
            settings.io.imageFile = [settings.io.imagePath, filesep,settings.io.imageName];
        else
            settings.io.imageName = char(unames(ii).fnames); % There is a mess how file names a read. Use char here for a .dat or .raw file
            settings.io.imageFile = char([settings.io.imagePath, filesep,settings.io.imageName]);
        end
    end
    
    settings.frc.note = [settings.io.imageName,'_bcg',num2str(settings.frc.bcgsub),'_sofis6lin'];
    
    if ~isempty(settings.io.roisx) && ~isempty(settings.io.roisy)
    settings.io.roit = {settings.io.roisx{ii},settings.io.roisy{ii}}; 
    end
    [sofi,sofi_c,settings,sofi_lin,stats,results] = SOFI2D_process(settings, str);
    
    if settings.io.matsave ==1
        save([settings.io.outputpath,filesep,settings.io.imageName,'.mat'],'sofi_c','sofi_lin');
%     save([settings.io.outputpath,filesep,settings.io.imageName,'.mat'],'sofi','sofi_c','sofi_lin','stats')
%      save([settings.io.outputpath,filesep,settings.io.imageName,'snr.mat'],'stats')
    end
    
    %Save SNR estimation results to separate file
    if settings.sys.jk ==1
    save([settings.io.outputpath,filesep,settings.io.imageName,'Jack_knife_SNR.mat'],'stats'); % write to a mat file
    end
    
end

c2=now;
set(0, 'DefaultFigureVisible', 'on');
disp(['Total time [hours]:',num2str((c2-c1)*24)])

set(0, 'DefaultFigureVisible', 'on');
%% Save settings and processing info
settings.io.proctime = (c2-c1)*24;

save([settings.io.outputpath,filesep,'config.mat'],'settings'); % write to a mat file

descfile_sofi3D(settings); % write to a text file
