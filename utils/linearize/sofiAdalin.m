function sofi_lin = sofiAdalin(sofi_c, stack, settings) 

wf = mean(sofi_c{1},3); % an average image i.e. equivalent to a widefield image
[params] = getimparams(wf);

wfbcg= max(0,wf-params.bcg); % subtract background
num = 10; % number of points used for the analysis

% find optimal linearization coefficients
[gcor,~,~,~] = kest_v6(sofi_c,num,wfbcg,stack(:,:,1:size(stack,3)));
% [gcor,rho_estB,coefs,parmap] = kest_v6b(sofi_c,num,wfbcg,stack(:,:,1:2000),settings);

gammas = mean(gcor,2,'omitnan')./(1:length(settings.sys.orders))';
gammas(1,:)=1;

% plot and save meantrace figure
if settings.io.figs ==1
    linada = mean(gammas,2);
    linst = 1./settings.sys.orders;
    figure,
    plot(linada,'bo--');hold on;
    text(1:length(linada),linada,num2str(round(1000*linada)./1000),'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
    plot(linst,'rx--');
    text(1:length(linst'),linst,num2str(round(1000*linst')./1000),'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','top')
    ylabel('Roots for SOFI linearization');
    xlabel('Cumulant order');
    xlim([1 length(linada)+1]);
    set(gca,'Xtick',1:length(linada))
end
% sofim=cellfun(@(x)squeeze(mean(reshape(x,size(x,1),size(x,2),5,[]),4)),sofi,'UniformOutput',0);
sofi_lin=sofiLinearizeCooseDeconvolution(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,settings, gammas);
% [sofi_lin]=sofiLinearize4(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,gammas,settings);