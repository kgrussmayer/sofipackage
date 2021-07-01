function sofi_lin = sofiAdalin(sofi_c, stack, settings) 

wf = mean(sofi_c{1},3); % an average image i.e. equivalent to a widefield image
[params] = getimparams(wf);

wfbcg= max(0,wf-params.bcg); % subtract background
num = 10; % number of points used for the analysis

% find optimal linearization coefficients
[gcor,~,~,~] = estimate_gamma_correction(...
    sofi_c,num,wfbcg,stack(:,:,1:size(stack,3)));

gammas = mean(gcor,2,'omitnan')./(1:length(settings.sys.orders))';
gammas(1,:)=1;

% plot and save meantrace figure
if settings.io.figs ==1
    linada = mean(gammas,2);
    linst = 1./settings.sys.orders;
    figure,
    plot(linada,'bo--');hold on;
    text(1:length(linada),linada,num2str(round(1000*linada)./1000),...
    'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
    plot(linst,'rx--');
    text(1:length(linst'),linst,num2str(round(1000*linst')./1000),...
    'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','top')
    ylabel('Roots for SOFI linearization');
    xlabel('Cumulant order');
    xlim([1 length(linada)+1]);
    set(gca,'Xtick',1:length(linada))
end

[sofi_lin]=sofiLinearize(sofi_c,settings.dec.fwhm,...
    settings.sys.orders,settings.dec.iter,settings, gammas);
% eof