function [ratio,density,dens,nummol] = molparams(sofi_c,sofi_lin,settings,shiftd)


%     if (~exist(outputPath,'dir'))
%     mkdir(outputPath);
%     end
%         [ratio,density,brightness]=sofiParametersMy(sofi_c);
        [ratio,density,brightness]=sofiParametersMy2(sofi_c);
%         [ratio,density,brightness]=sofiParameters(sofi_c);
        
        sofimask = mean(sofi_lin{4},3);
        alpha = sofimask;
        alpha_max = max(alpha(:));
        alpha = imadjust(alpha./alpha_max,[settings.molpar.thresh 1],[0 1]);

        % alpha = alpha./max(alpha(:));
        alpha(alpha>0) = 1;
        density(isinf(density)) = 0;

        I_ratio = ratio.*sofimask;
        I_bright = brightness.*alpha;
        I_density = density.*sofimask;

        bsofi=sofiBalance(sofi_lin,ratio);
           
% save([settings.io.outputpath,filesep,settings.io.imageName,'_molparams.mat'],'ratio','density','brightness','sofimask','bsofi',...
%             'I_ratio','I_bright','I_density'); 
% saveImageStack(bsofi,settings.io.outputpath,[settings.io.imageName,'bsofi'],[],16)

%% number of molecules, depends on max and min threshold - too low values affected by low SNR, to high values are outliers
% roi{2} = density(436:436+80,68:68+80);
% roi{3} = density(316:316+80,244:244+80);
% roi{4} = density(196:196+80,432:432+80);
% xaxis = linspace(0.05,0.1,5);
xaxis = [0, 0.01, 0.1, 0.5, 1, 1.5, 2];

clear dens
clear nummol;
% density = medfilt2(density,[2 2]);
count = 1;
% for th = xaxis
%     i = find(density>th & density<300);
%     i = find(density>th);
%     i = find(density>th & I_density<200);
%     dens(count) = mean(density(i)).*(1/(settings.frc.pixelsize/1000/4).^2);
%     figure, imshow(I_density,[]);colormap(jet)
%     figure, imshow(bsofi,[]);colormap(jet);hold on;
    rx = 90; %+ shiftd; %37
    ry = 90;  %%roi2 88, rlen 58
    rlen = 40; %%roi3 111, 37, roi4 100, 40 
    % roi5 98 37
    % roi6 95 40 
    % roi7 90 40 
    %roi8 88 37
    roix = rx+1:rx+rlen;
    roiy = ry+1:ry+rlen;
%     plot(roix,repmat(ry+1,length(roix),1),'k');
%     plot(roix,repmat(ry+rlen,length(roix),1),'k');
%     plot(repmat(rx+1,length(roix),1),roiy,'k');
%     plot(repmat(rx+rlen,length(roix),1),roiy,'k');
    
%     mid = floor(size(density)./2);
    density_roi = medfilt2(density(roiy,roix),[2 2]);
%     densnoise = density;
%     densnoise(roiy,roix) = NaN;
%     densnoise = std(densnoise(:), 'omitnan');
%     densnoise = median(densnoise(:));
%     density_roi = max(density_roi-densnoise,0);
%     density_roi(density_roi==0)=NaN;
%     i = find(density_roi>densnoise);
    dens(count) = mean(density_roi(:)).*(1/(settings.frc.pixelsize/1000/4).^2);
    
%     density_roi = density(mid(1)-5*4:mid(1)+4*4,mid(2)-5*5:mid(2)+4*4);
%     density_roi = density(35:35+40,35:35+40);
%       density_roi = density(30:30+40,30:30+40);
%     dens(count) = mean(density_roi(:)).*(1/(settings.frc.pixelsize/1000/4).^2); 
%     nummol(count) = sum(density(i)); 
%     i = find(density_roi>1);
    nummol(count) = sum(I_density(:)); 
    count = count +1;
% end

%% show figures
if settings.io.figs ==1   
    figure, 
    imshow(I_ratio,[0 0.3]);colormap('jet'); title('On-time ratio');colorbar;
    if settings.io.figsave ==1
        save2pdf([settings.io.outputpath,filesep,settings.io.imageName,'_onratio.pdf'],gcf,300);
        saveFigure(gcf,settings.io.outputpath,[settings.io.imageName,'_onratio'],'png');
    end

    I_dens = I_density;
    % I_dens(I_dens>10)=10;

    figure, 
    im = nperc(medfilt2(I_density,[3 3]),0.1);
    imshow(im,[]);colormap('jet'); title('Density map');colorbar;

    ylabel('Molecular density [molecules/px^2]')
    set(gca,'yaxislocation','right');
    if settings.io.figsave ==1
        save2pdf([settings.io.outputpath,filesep,settings.io.imageName,'_density_v2.pdf'],gcf,300);
        saveFigure(gcf,settings.io.outputpath,[settings.io.imageName,'_density_v2'],'png');
    end

    figure,
    I_density_au = medfilt2(I_density)./max(max(medfilt2(I_density)));

    im = nperc(I_density_au,0.1);
    imshow(im,[]);colormap('jet'); title('Density map');colorbar;
    ylabel('Molecular density [a.u.]')
    set(gca,'yaxislocation','right');
    if settings.io.figsave ==1
        save2pdf([settings.io.outputpath,filesep,settings.io.imageName,'_density_v1.pdf'],gcf,300);
        saveFigure(gcf,settings.io.outputpath,[settings.io.imageName,'_density_v1'],'png');
    end

end

%     figure, 
%     subplot(311);imshow(I_ratio,[]);colormap('jet'); title('On time ratio');colorbar;
%     subplot(312);imshow(I_bright,[]);colormap('jet'); title('Brightness');colorbar;
%     subplot(313);imshow(I_density,[]);colormap('jet'); title('Number of molecules'); colorbar;
