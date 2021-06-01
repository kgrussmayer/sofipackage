function results = estimResol(sofi_lin, settings)
% estimate resolution

%% FRC

%
%%% calculate FRC
%
count = 1;
pxsize =  settings.frc.pixelsize;

for order = settings.frc.orders
  
    [im1, im2] = prepareFrc(sofi_lin,settings,order);

    dip_initialise; % initialize dip image package

%     if roi ==1
%     % crop region of interest [roiy, roix]
%     im1 = im1(settings.roix ,settings.roix);
%     im2 = im2(settings.roix ,settings.roix);
%     end
    settings.frc.pixelsize = pxsize./order;
    [resMid, resHigh, resLow, frcCurves] = frcs(im1,im2, settings.frc);

    results(count).resHigh = resHigh;
    results(count).resMid = resMid;
    results(count).resLow = resLow;
    results(count).frc_curves = frcCurves;
    
    %
    %%% plot FRC
    %

    if settings.io.figs ==1   
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

        saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
            [settings.io.imageName,'_frc_sofi',num2str(order)],settings.io.figformat,0);
    end
    
    count = count+1;
end
settings.frc.pixelsize = pxsize;
% save([settings.io.outputpath,filesep,'FRC_bcgsub_',num2str(settings.frc.bcgsub),'.mat'],'results','settings','-v7.3');