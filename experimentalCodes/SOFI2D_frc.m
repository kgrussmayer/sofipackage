% 2D SOFI processing: calculate FRC from precomputed data
% Tomas Lukes, tomas.lukes@epfl.ch 

config;

% settings.sys.orders = 4;
% settings.frc.orders = settings.sys.orders;

settings.io.figs =0;
settings.io.outputpath = [imagePath,filesep,'postproc',filesep,fnames{1},'_wsize',num2str(settings.sys.wsize),'_fwhm',num2str(settings.dec.fwhm),'bl',...
    num2str(settings.io.blcor ),'dc',num2str(settings.io.dcor),'_start',num2str(settings.sys.start)];% output folder for results

settings.io.imageName = fnames{1}; 
settings.io.imageFile = [imagePath, filesep,settings.io.imageName];
    
%
%%% Load mat file
%

% [stack,frames] = loadStack(settings);
% settings.sys.frames = frames;
% numim = size(stack,3);
% maxf = 4; %first 1000 frames
% c{settings.sys.orders}=double(stack(:,:,:)); 
load([settings.io.imageFile,'.mat']);

%% Deconv
% settings.dec.fwhm = 2;
% settings.dec.iter = 5;
% 
% sofi=sofiLinearize(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter);

%%
%
%%% FRC estimation
%

settings.frc.run = 1;
settings.frc.orders = 3;
settings.frc.bcgsub = 1;%1.3;
settings.frc.pixelsize = settings.sys.pxy;

if settings.frc.run == 1
    results = estimResol(sofi, settings);
end

%%
figure,
for ii = 1:7
    subplot(1,7,ii);
    plot(results.frc_curves{1, ii})
    hold on;
    plot(1:282,1/7*ones(1,282),'r--')
end




resMid = results(1).resMid;
figure, 
        plot(1:length(resMid),resMid,'x--');
        xlabel('Sector number');
        ylabel('sFRC [nm]');
        % set(gca,'XTick',[1:7])
        xframes = {[1:length(resMid)-1],'full FRC'};
        set(gca,'XTickLabel',xframes);

        f1 = @(x) sprintf('   %.1f',x); 
        str = cellfun(f1, num2cell(resMid), 'UniformOutput', false);
        text(1:length(resMid),resMid,str,'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','bottom')