%% This file reproduce the behavior of SOFI2D_batch
% in a cell structure. This allows for easy test prior to batch processing
% and custom for loop over specific part of the sofi processing
% 

addpath(genpath('utils'));
config;

rootPath = 'C:\Users\adrien\Desktop\sofi_blinking';
cDir = cd; cd(rootPath);
[fname,pname] = uigetfile({'*.*','All files'});
cd(cDir)

disp(['Path name : ',pname])
disp(['File name : ',fname])

ind = strfind(fname,'.');

% define settings structure
settings.io.fileExt = fname(ind(end)+1:end);
settings.io.imagePath = pname;
settings.io.fname = fname;
settings.io.imageName = fname(1:ind(end)-1);

%% Processing settings
% bleaching correction
settings.io.blcor = 0;

% drift correction
settings.io.dcor = 0;
settings.dcor.type = 'SOFI';
settings.dcor.drift = [];

% typical variables for batch processing
settings.sys.wsize = 451;
settings.sys.start = 1;
settings.sys.end = [];
settings.sys.frames = [];
settings.sys.avg = 1; % average sofi subsequence after drift correction for faster linearization

% linearization
settings.dec.lin = 1; % set to 1 if you want to linearize SOFI
settings.dec.fwhm = 2.5;
settings.dec.iter = 10;
settings.dec.lincoef = 0.7;

% deconvolution
settings.augLag.deconv = 0; % set to 1 if you want to deconvolve the linerized SOFI
settings.augLag.gamma = 2000;
settings.augLag.FWHM = 3;
settings.augLag.maxIter = 20;

%% data Loading
data = load_tifFile([settings.io.imagePath, filesep, settings.io.fname]);
settings.sys.frames = size(data,3);

%% looping parameter
field = 'wsize'; %'wsize','start','end' ; 'fwhm','iter','lincoef' ; 'gamma','FWHM','maxIter' ; 'blcor','dcor'
param = {451}; % {50 200 400 600 800 1000 2000 4000 []};
% param = {1 200 400 600 800 1000 1200 1400 1600 1800 2000};
settings.io.outputpath = [settings.io.imagePath,filesep,'results_',getID(5),'_field_',field,'_',settings.io.imageName];

% store input settings
s = settings;

sofi_cf = cell(4,1);    sofi_decf = cell(4,1);
sofi_f = cell(4,1);     sofi_linf = cell(4,1);

%%
for k = 1:length(param)
    disp(['---- field : ',field,' ----'])
    disp(['Field value : ',num2str(param{k})])
    
    settings = s;
    if isempty(settings.sys.frames); settings.sys.frames = size(data,3); end
    
    % set settings according to the specified field 
    settings = setnestedfield(settings,field,param{k});
    
    stack = data;
    
        % crop stack 
	if isfield(settings.io,'roi')
        if ~isempty(settings.io.roi)
            roi = settings.io.roi;
            stack = stack(roi{2},roi{1},:);
        end
    end
    % check and apply settings
    if ~isempty(settings.sys.end); stack = stack(:,:,1:settings.sys.end); end
    if ~isempty(settings.sys.start); stack = stack(:,:,settings.sys.start:end); end
    if isempty(settings.sys.wsize); settings.sys.wsize = size(stack,3); end
    
% bleaching correction
    if settings.io.blcor == 1
        [settings,results] = plotmitrace(stack,settings);
        [stack,results] = bleachcor(stack, settings,results);
    end

    % Ton estimate
    if settings.ton.run ==1
        results = estimateton(stack,settings,results);
    end

% Drift correction 
    % drift correction based on tracking fiducial markers
    if settings.io.dcor == 1 
        if strcmp(settings.dcor.type,'TS')
            stack = driftcorTS(stack, settings,frames);
        elseif strcmp(settings.dcor.type,'LBEN_PALM')
            stack = driftcor(stack,settings);
        elseif strcmp(settings.dcor.type,'SOFI')
            disp('Drift correction from stack')
        else
            disp(['Wrong drift corr type :,',settings.dcor.type])
        end
    end
    
    % SOFI cumulants
    disp('Cumulants')
    [sofi_c, settings,stats] = sofiCumulantsAll(stack, settings);
    

    % post-processing drift correction
    if settings.io.dcor == 1 && strcmp(settings.dcor.type,'SOFI')
        im_fix = sofi_c{2}(:,:,1); 
        
        [sy,sx,~] = size(im_fix);
        [xx, yy]=meshgrid(1:sx,1:sy);
        
        disp('Drift estimation')
        disp(['|',repmat('-',[size(sofi_c{2},3) 1])','|'])
        fprintf('%s','|')
        for ii = 1:size(sofi_c{2},3)
            
            im_mov = sofi_c{2}(:,:,ii);
            [dx, dy] = ccrShiftEstimation(im_fix,im_mov,10);
            results.sofidrift(:,ii) = [dx, dy]./2; % convert the shift in
            
%             if ii > 1 && settings.io.dcor == 1 
%                 % drift correction for 2nd order SOFI images
%                 imReg=interp2(xx,yy,sofi_c{2}(:,:,ii),xx+dx,yy+dy,'linear');
%                 imReg(isnan(imReg)) = 0;
%                 sofi_c{2}(:,:,ii) = imReg;
%             end
            fprintf('%s','-')
        end
        fprintf('%s\n','|')
         settings.dcor.drift = results.sofidrift;
         
         figure(623423)
         plot(-(settings.sys.pxy/2)*results.sofidrift(1,:),'x-','linewidth',2);hold on;
         plot(-(settings.sys.pxy/2)*results.sofidrift(2,:),'+-','linewidth',2)
         xlabel('SOFI2 subwindows'); ylabel('Drift [WF pixels]');legend('X drift','Y drift'); hold off
         saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],'SOFI2_drift_estimate','.png')
         
        % apply drift correction
     	 for n = 1:length(sofi_c)
            s = sofi_c{n};
            for h = 1:size(sofi_c{n},3)
                t0 =(h-1)*settings.sys.wsize+1;
                t1 = (h)*settings.sys.wsize;
                if n == 1
                    sofi_c{n}(:,:,h) = imtranslate(mean(stack(2:end-2,2:end-2,t0:t1),3),settings.dcor.drift(:,h)');
                else
                    sofi_c{n}(:,:,h) = imtranslate(sofi_c{n}(:,:,h),n.*settings.dcor.drift(:,h)');
                end
            end
        end
    end


    % average the cumulants for faster post-processing
    if settings.sys.avg == 1
    disp('Averaging raw cumulants')
        temp = cell(4,1);
        for n = settings.sys.orders
            temp{n} = mean(sofi_c{n},3);
        end
        sofi_c = temp;
    end    
    if settings.dec.lin
        % Adaptive linearization
        if max(settings.sys.orders) >= 4
            disp('Adaptive linearization')
            sofi_lin = sofiAdalin(sofi_c, stack, settings);
        end
        % Linearization
        disp('Regular linearization')
        sofi=sofiLinearize(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,settings);
    
        % Deconvolution
        sofi_dec = [];
        if settings.augLag.deconv
        disp('Augmented lagrangien deconvolution')
            for n = settings.sys.orders
                [sy,sx,Nw] = size(sofi{n});
                h2 = fspecial('Gaussian',[sy,sx],settings.augLag.FWHM*sqrt(n)/2/sqrt(2*log(2)));
                for h = 1:Nw
                    tmp=(sofi{n}(:,:,h));
                    tmp(tmp<0)=0;
                    scale = max(tmp(:));
                    sofi_dec{n}(:,:,h) = scale*deconvAugLag3D(tmp/max(tmp(:)),fftn(h2,[sy,sx]),settings.augLag);
                end
            end
        end
    end
    % storing results
    for n = settings.sys.orders
        sofi_cf{n}(:,:,k) = mean(sofi_c{n},3);
        if settings.dec.lin
            sofi_linf{n}(:,:,k) = mean(sofi_lin{n},3); 
            sofi_f{n}(:,:,k) = mean(sofi{n},3);
            if settings.augLag.deconv 
                sofi_decf{n}(:,:,k) = mean(sofi_dec{n},3);
            end
        end
    end
end

disp('End of processing')

%% Saving results
disp('Saving results')
out = mkdir(settings.io.outputpath);
settings = setnestedfield(settings,field,param);
settings.field = field;
settings.fieldValue = param;
saveSettingsTxt(settings,[settings.io.outputpath,filesep,'000_settings.txt'])
    
outputpath = settings.io.outputpath;
fname = settings.io.fname;
ind = strfind(fname,'.'); fname= fname(1:ind(end)-1);
for n = settings.sys.orders
    writeTIFF(sofi_cf{n},[outputpath,filesep,'sofi_c',num2str(n),'_',fname])
    if settings.dec.lin
        writeTIFF(sofi_linf{n},[outputpath,filesep,'sofi_lin',num2str(n),'_',fname])
        writeTIFF(sofi_f{n},[outputpath,filesep,'sofi_',num2str(n),'_',fname]) 
        if settings.augLag.deconv 
            writeTIFF(sofi_decf{n},[outputpath,filesep,'sofi_dec',num2str(n),'_',fname])
        end
    end
end

% % save SOFI stacks (sofi_xxxf are averaged)
% outputpath = 'L:\segaud\simulMTF\simulations\results\1MTF_3rep_0o60o120o_3phases_500frames_sofisim\216nmgrating_SOFI_stacks';
% fname = settings.io.fname;
% ind = strfind(fname,'.'); fname= fname(1:ind(end)-1);
% 
% for n = settings.sys.orders
%     writeTIFF(sofi_c{n},[outputpath,filesep,'sofi_c',num2str(n),'_',fname]) 
%     writeTIFF(sofi_lin{n},[outputpath,filesep,'sofi_lin',num2str(n),'_',fname])
%     writeTIFF(sofi{n},[outputpath,filesep,'sofi_',num2str(n),'_',fname]) 
%     writeTIFF(sofi_dec{n},[outputpath,filesep,'sofi_dec',num2str(n),'_',fname]) 
% end
% for n = settings.sys.orders
%     tempSofi=[];
%     tempSofi(:,:,1)=mean(sofi_c{n}(:,:,1:2:end),3);
%     tempSofi(:,:,2)=mean(sofi_c{n}(:,:,2:2:end),3);
%     writeTIFF(tempSofi,[outputpath,filesep,'sofi_c',num2str(n),'_',fname])
%     tempSofi=[];
%     tempSofi(:,:,1)=mean(sofi_lin{n}(:,:,1:2:end),3);
%     tempSofi(:,:,2)=mean(sofi_lin{n}(:,:,2:2:end),3);
%     writeTIFF(tempSofi,[outputpath,filesep,'sofi_lin',num2str(n),'_',fname])
%     tempSofi=[];
%     tempSofi(:,:,1)=mean(sofi{n}(:,:,1:2:end),3);
%     tempSofi(:,:,2)=mean(sofi{n}(:,:,2:2:end),3);
%     writeTIFF(tempSofi,[outputpath,filesep,'sofi_',num2str(n),'_',fname]) 
%     tempSofi=[];
%     tempSofi(:,:,1)=mean(sofi_dec{n}(:,:,1:2:end),3);
%     tempSofi(:,:,2)=mean(sofi_dec{n}(:,:,2:2:end),3);
%     writeTIFF(tempSofi,[outputpath,filesep,'sofi_dec',num2str(n),'_',fname]) 
% end

disp('EOF')