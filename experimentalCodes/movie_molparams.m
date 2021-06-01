clear all; close all;

config;

fp = 'C:\Users\lukes\Documents\MEF_10ms532nm_FP_ND2_405_50mw_ND1_DU897_BV_2576_wsize500_fwhm2.5bl0dc1_start1';
fn = 'MEF_10ms532nm_FP_ND2_405_50mw_ND1_DU897_BV_2576';
load([fp,filesep,fn]);

% Calculate molecular parameters
settings.io.outputpath = fp; 
settings.io.imageName = fn;
settings.io.figsave = 0;
settings.molpar.thresh = 0.2;
settings.io.figsave =1;


%%

frames = 4000;
fps = 100; 
wsize = frames/settings.sys.wsize;

vx = 150;
vy = 90;

avgshift_x = settings.sys.wsize*(vx / (60*fps)); % average shift in nm
avgshift_x = avgshift_x ./(settings.sys.pxy./4); % avg shift in px of sofi4

avgshift_y = settings.sys.wsize*(vy / (60*fps)); % average shift in nm
avgshift_y = avgshift_y ./(settings.sys.pxy./4); % avg shift in px of sofi4

im = mean(sofi{4},3);

[im_roi, roi] = imcrop(im,[]); % top left corner x, top left corner y, width x, width y
roix = roi(1):roi(1) + roi(3);
roiy = roi(2):roi(2) + roi(4);

%%
for ii = 1:size(sofi_c{2},3)-wsize+1;
  
sofi_cw{1} = sofi_c{1}(:,:,ii:ii+wsize-1);    
sofi_cw{2} = sofi_c{2}(:,:,ii:ii+wsize-1);
sofi_cw{3} = sofi_c{3}(:,:,ii:ii+wsize-1);
sofi_cw{4} = sofi_c{4}(:,:,ii:ii+wsize-1);

roix = round(roix + (ii-1)*avgshift_x);
roiy = round(roiy + (ii-1)*avgshift_y);

[ratio,density,dens,nummol] = molparams2roi(sofi_cw,sofi_lin,settings,roix,roiy,ii);
out.dens(ii) = dens;
% out.nummol(ii) = nummol;
out.roiy{ii} = roiy;
out.roix{ii} = roix;
dmaps(:,:,ii) = density;
end

%%
figure, bar(out.dens);
xlabel('Frame number');
ylabel('Density in the ROI');

save2pdf([settings.io.outputpath,filesep,settings.io.imageName,'_densityMovie','.pdf'],gcf,300);
saveFigure(gcf,settings.io.outputpath,[settings.io.imageName,'_densityMovie'],'png');
        
save([settings.io.outputpath,filesep,settings.io.imageName,'_livecell_rois_dens','.mat'],'out','dmaps');
%%
figure, 
imshow(im,[]);hold on;
plotRect(round(roi),'ROI','red')

%%
figure, 
imshow(dmaps(:,:,1),[0 5]);colormap('jet');

