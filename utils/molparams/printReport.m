function printReport(settings,results,sofi_c,sofi_lin,imIn,mparams)
titleFontSize = 8;

f1 = figure('units','normalized','outerposition',[0 -0.08 0.9 1.2]);

% Widefield
s1 = subplot(8,3,[1 4]);
imshow(imIn./max(imIn(:)),[0 1]);
title('Widefield','FontSize', titleFontSize);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;

% Meantrace

s2 = subplot(8,3,[2 5]);
plot(results.blcor.xaxis',results.blcor.mitrace);
hold on;
plot(results.blcor.xaxis,results.blcor.fitcurve,'k');
plot(results.blcor.xaxis,results.blcor.mitrace_cor);
ylabel('Mean intensity');
title(['Est. bleaching lifetime ',num2str(results.blcor.fitinfo.b),' frames'],...
    'FontSize', titleFontSize);
ylim([0.3 1.2]);
legend('meantrace','exp fit','bleach. corrected','Location','southwest');
set(gca,'FontSize',8);

% Autocorrelation
s3 = subplot(8,3,[3 6]);
plot(1:length(results.blcor.corelf),results.blcor.corelf);hold on;
%     plot(1:length(results.blcor.corela_cor),results.blcor.corela_cor);
ylabel('Autocorrelation'); legend('input data','bleach. corrected');
set(gca,'FontSize',8);
freezeColors;

% On-time ratio map
s4 = subplot(8,3,[7 10]);
imshow(mparams.I_ratio,[0 0.5]);
psub = get(s4, 'position');
colormap('jet');colorbar('Position',[psub(1)+0.73*psub(3), psub(2), psub(3)/20, psub(4)]);
title('On-time ratio','FontSize', titleFontSize);
set(gca,'FontSize',8);
freezeColors;

% Density map
s5 =subplot(8,3,[8 11]);
I_density_au = medfilt2(mparams.I_density)./max(max(medfilt2(mparams.I_density)));
im = nperc(I_density_au,0.1);
imshow(im,[]);
psub = get(s5, 'position');
colormap('jet');colorbar('Position',[psub(1)+1.05*psub(3), psub(2), psub(3)/20, psub(4)]);
title('Density map','FontSize', titleFontSize);
%     ylabel('Molecular density [a.u.]')
%     set(gca,'yaxislocation','right');
set(gca,'FontSize',8);
freezeColors;

% On-time ratio histogram
subplot(8,3,[9]);
temp = mparams.ratio(:);
temp(temp==0)=[];
hist(temp,[0:0.01:1],'b');
title('On-time ratio histogram','FontSize', titleFontSize);
xlim([0 1]);colormap('jet');
set(gca,'FontSize',8);
freezeColors;

% Density histogram
subplot(8,3,[12]);
temp = I_density_au(:);
temp(temp==0)=[];
hist(temp,[0:0.01:1]);title('Normalized density histogram','FontSize', titleFontSize);
xlim([0 1]);colormap('jet');
set(gca,'FontSize',8);
freezeColors;

% raw SOFI
subplot(8,3,[13 16]);imshow(mean(sofi_c{2},3),[]),title('Raw SOFI order: 2','FontSize', titleFontSize);
subplot(8,3,[14 17]);imshow(mean(sofi_c{3},3),[]),title('Raw SOFI order: 3','FontSize', titleFontSize);
subplot(8,3,[15 18]);imshow(mean(sofi_c{4},3),[]),title('Raw SOFI order: 4','FontSize', titleFontSize);

% bSOFI
subplot(8,3,[19 22]);imshow(mean(sofi_lin{2},3),[]),title('Lin SOFI order: 2','FontSize', titleFontSize);
subplot(8,3,[20 23]);imshow(mean(sofi_lin{3},3),[]),title('Lin SOFI order: 3','FontSize', titleFontSize);
subplot(8,3,[21 24]);imshow(mean(sofi_lin{4},3),[]),title('Lin SOFI order: 4','FontSize', titleFontSize);
colormap('morgenstemning');freezeColors;

colormap(s4,'jet');
colormap(s5,'jet');

tightfig;

if settings.io.figsave ==1
    fnameout = [settings.io.outputpath,filesep,'figs',filesep,settings.io.imageName,'_report'];
    save2pdf([fnameout,'.pdf'],gcf,150);
end
% eof