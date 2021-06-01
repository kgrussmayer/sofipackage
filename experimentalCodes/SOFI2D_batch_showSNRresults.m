% 2D SOFI - batch process all input files
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2015a

% Load file names and processing settings from the config file
config;
% config_local;

c1=now;
warning('off','all')
addpath(genpath('utils'));

if settings.io.figshow ==1
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');  
end
    
%% SOFI 2D batch 
[fnames2,~] = getFileNums(fnames); % detect numbers from input image names 
unames = uniNames(fnames,fnames2); % group unique names

for ii = 1:numel(fnames)
    c(ii)=now;
    disp(['Processing file number: ',num2str(ii),' from ',num2str(numel(fnames))]);
    
    if settings.io.concatOn == 1;
        settings.io.imageName = unames(ii).fnames;  
    else
        settings.io.imageName = fnames{ii}; 
        settings.io.imageFile = [settings.io.imagePath, filesep,settings.io.imageName];
    end
    
    settings.frc.note = [settings.io.imageName,'_bcg',num2str(settings.frc.bcgsub),'_sofis6lin'];
    
    if ~isempty(settings.io.roisx) && ~isempty(settings.io.roisy)
    settings.io.roit = {settings.io.roisx{ii},settings.io.roisy{ii}}; 
    end
    
    [sofi,sofi_c,settings,sofi_lin,stats,results] = SOFI2D_process(settings);
    
    if settings.io.matsave ==1
        save([settings.io.outputpath,filesep,settings.io.imageName,'.mat'],'sofi_c','sofi_lin');
%     save([settings.io.outputpath,filesep,settings.io.imageName,'.mat'],'sofi','sofi_c','sofi_lin','stats')
%      save([settings.io.outputpath,filesep,settings.io.imageName,'snr.mat'],'stats')
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


%% Show SNR results

if settings.sys.jk == 1 % show the SNR only if the calculation was performed
    
    if isstruct(stats)
    snrs = stats.s;
    vars = stats.v;
    bias = stats.b;   

%     snrs=sofiAllFlatten(snrs,settings.sys.orders);

    for ii = 1:size(snrs{2},3);
    % snrsm(:,:,ii) = medfilt2(snrs{2}(:,:,ii),[3 3]);
    snrsm(:,:,ii) = imfilter(snrs{2}(:,:,ii),fspecial('gaussian',[5 5],1));
    end

    blockSize = 2;

    vars2 = binStack(abs(vars{2}), blockSize);
    snrs2 = binStack(abs(sofi_c{2}), blockSize)./sqrt(vars2);
%     snrs2 = abs(sofi_c{2})./sqrt(abs(vars{2}));
        figure, 
        subplot(131);
        imshow(imresize(mean(snrs{2}(1:2:end,1:2:end,:),3),2),[]);colormap('jet');colorbar;
        subplot(132);
        imshow(mean(snrs2(:,:,:),3),[]);colormap('jet');colorbar;
        subplot(133);
        imshow(mean(snrsm(:,:,:),3),[]);colormap('jet');colorbar;
    
    end

    % snrs2 = sofi_c{2}./sqrt(vars{2});
    % snrs2(isnan(snrs2)) = 0;
    snr = squeeze( mean(mean(snrsm,1,'omitnan'),2,'omitnan'));
    snr2 = squeeze( mean(mean(snrs2,1,'omitnan'),2,'omitnan'));
    
    figure, 
    plot(snr);hold on;
    plot(snr2);
    
    %
    %%% SNR histogram
    %
    
    figure('position',[300 300 1200 300]);
%     for ii = 1:size(snrs{2},3)

    % temp = mean(snrs2(:,:,1),3);
%     temp = mean(snrs2(:,:,:),3);
% temp = imresize(mean(snrs{2}(1:2:end,1:2:end,:),3),2);
    ii = 1;
    temp = mean(snrs{2}(1:2:end,1:2:end,:),3);
    temp(isinf(temp))=0;

    subplot(131);imshow(temp,[]);colormap('jet')

    disp(mean(temp(:),'omitnan'));
    SNR_a(ii) = mean(temp(:),'omitnan'); % average SNR
    temp = round(temp*100)./100;

    subplot(132);
    hist(temp(:),length(unique(temp)));
    [counts,centers] = hist(temp(:),length(unique(temp)));hold on;
    % plot(counts);
    % hist(round(temp(:)),length(unique(round(temp))));


    %         mitrace = mitrace./max(mitrace);
            s = fitoptions('Method','NonlinearLeastSquares',...
                           'Lower',[0,0,0],...
                           'Upper',[Inf,Inf,Inf],...
                           'Startpoint',[1 1 1]);
            f = fittype('a*exp(-x/b)+c','options',s);

            [c2,gof2] = fit((0:length(counts)-1)',counts(1:end)','gauss2');
    %         ft=fit(x,y,'gauss2') 
            a1 = c2.a1;
            b1 = c2.b1;
            c1 = c2.c1;

            a2 = c2.a2;
            b2 = c2.b2;
            c2 = c2.c2;

            SNR_b(ii) = centers(round(b1)); % average bacground SNR
            SNR_f(ii) = centers(round(b2)); % average foreground SNR
            thresh = graythresh(temp./max(temp(:))).*max(temp(:));
            [minv, minid] = min(abs(centers-thresh));
            plot(ones(1,max(counts)).*thresh,1:max(counts),'r--');
            ylabel('SNR histogram');xlabel('SNR');
            title(['SNR threshold = ',num2str(thresh)]);
            
            subplot(133);
            
            x=(0:length(counts)-1);
            plot(x',counts);hold on;
            plot(a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2),'k');
            ylabel('SNR histogram');xlabel('Bin number');
%             plot(ones(1,max(counts)).*minid,1:max(counts),'r--');
%             set(gca,'XTick',x(1:100:end))
%             xaxis = round((0:max(temp(:))/(length(counts)-1):max(temp(:))).*10)./10;
%             set(gca,'XTickLabel',[xaxis(1:100:end)])
            hold off;
            title(['SNR_f =',num2str(SNR_f(ii)),'  SNR_b =',num2str(SNR_b(ii))]); 
    %         saveFigure(gcf,settings.io.outputpath,[settings.io.imageName,'SNR_fit_f',num2str(ii)],'png',0)
%     end
%%
       
      mask = im2bw(temp./max(temp(:)),graythresh(temp./max(temp(:))));
      figure,
%       mask = zeros(size(temp));
%       mask(temp>4)=1;
      
      imshow(mask,[]);
%     figure, 
%     plot(SNR_a);hold on;
%     plot(SNR_b);
%     plot(SNR_f);
%     legend('avg SNR','avg SNR background','avg SNR foreground'); 
    % saveFigure(gcf,settings.io.outputpath,[settings.io.imageName,'SNRs'],'png',0)

end

