function [mparams] = molparams(sofi_c,sofi_lin,settings)
% wrapper function for sofi molecular parameters calculation

% estimate molecular parameters
[ratio,density,brightness]=sofiParameters(sofi_c);

% molecular parameters are relevant only for pixels with high signals
% create mask using linearised 4th order SOFI
sofimask = sofi_lin{4};
alpha = mean(sofimask,3);
alpha_max = max(alpha(:));
alpha = imadjust(alpha./alpha_max,[settings.molpar.thresh 1],[0 1]);

alpha(alpha>0) = 1;
density(isinf(density)) = 0;

% possibly create mask using balanced SOFI
bsofi=sofiBalance(sofi_lin,ratio);

% masking pixels
I_ratio = ratio.*alpha;
I_bright = brightness.*alpha;
I_density = density.*alpha;

mparams = makestruct(ratio,brightness,density,I_ratio,I_bright,...
    I_density,alpha,bsofi);

save([settings.io.outputpath,filesep,settings.io.imageName,...
    '_molparams.mat'],'ratio','density','brightness','sofimask',...
    'bsofi','I_ratio','I_bright','I_density');
% save results matlab structure to matfile
saveImageStack(bsofi,settings.io.outputpath,...
    [settings.io.imageName,'bsofi'],[],16);

% show figures
if settings.io.figs ==1
    figure,
    imshow(I_ratio,[0 0.5]);colormap('jet'); title('On-time ratio');
    colorbar;
    if settings.io.figsave ==1
        save2pdf([settings.io.outputpath,filesep,'figs',filesep,...
            settings.io.imageName,'_onratio.pdf'],gcf,300);
        saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
            [settings.io.imageName,'_onratio'],'png');
    end
    
    figure
    im = nperc(medfilt2(I_density,[3 3]),0.1);
    imshow(im,[]);colormap('jet'); title('Density map');colorbar;
    
    ylabel('Molecular density [molecules/px^2]')
    set(gca,'yaxislocation','right');
    if settings.io.figsave ==1
        save2pdf([settings.io.outputpath,filesep,'figs',filesep,...
            settings.io.imageName,'_density_v2.pdf'],gcf,300);
        saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
            [settings.io.imageName,'_density_v2'],'png');
    end
    
    figure
    I_density_au = medfilt2(I_density)./max(max(medfilt2(I_density)));
    
    im = nperc(I_density_au,0.1);
    imshow(im,[]);colormap('jet'); title('Density map');colorbar;
    ylabel('Molecular density [a.u.]')
    set(gca,'yaxislocation','right');
    if settings.io.figsave ==1
        save2pdf([settings.io.outputpath,filesep,'figs',filesep,...
            settings.io.imageName,'_density_v1.pdf'],gcf,300);
        saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
            [settings.io.imageName,'_density_v1'],'png');
    end
    
end
% eof