%%% Test Fourier ring correlation metric on SOFI images, calculate FRC in
%%% the same way as Henrik on Palm images, 

%%% Tomas Lukes
%%% created 23.9.2015
% 
% clear all;
% close all;

function frcs_03(sofi_c, sofi, imagePath, settings,pixelsize)
% im1 = logical(round(rand(1,round(size(sofi_lin{2},3)/2))));

% single tif file
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq500_fwhm3_mov11_roi_1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq500_fwhm8_20150619_fixed_cell_meos2_mov_11_roi_1_dc2x';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq500_fwhm4_20150619_fixed_cell_meos2_mov_11_roi_1_dc';

% imagePath = 'seq500_fwhm2_20150619_fixed_cell_meos2_mov_11_roi_1_dc';
% imagePath = 'seq500_fwhm2_20150619_fixed_cell_meos2_mov_11_roi_1_nodc';

% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq500_fwhm3_mov3b_roi_1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize500_fwhm4bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize250_fwhm4bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize500_fwhm5bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize500_fwhm6bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize1000_fwhm4bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize1000_fwhm5bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize2000_fwhm4bl1dc1';
% imagePath = 'H:\tlukes\data_Hendrik_livecells\seq_mov11_roi1_wsize5000_fwhm4bl1dc1';

% outputPath = 'H:\tlukes\data_Hendrik_livecells\seq500_fwhm3_mov11_roi_1\frc_measurements';
% outputPath = 'H:\tlukes\data_Hendrik_livecells\FRC_mEos_mov11\';
% outputPath = [imagePath,filesep,'FRC measurements 3d_sec'];
outputPath = [imagePath];
% 
if ~exist(outputPath)
    mkdir(outputPath)
end

% imageName = '20150619_fixed_cell_meos2_paxillin_fiducial_mov_11';
% imageName = '20150619_fixed_cell_meos2_paxillin_fiducial_mov_11sofi_lucy2Movie2';
% imageName = '20150722_mov_3b';

% path with folder which contains numbers of frames within each subset - in
% mat files from Hendrik
% frcPath = 'H:\tlukes\data_Hendrik_livecells\for_tomas1\fourier ring correlation\20150619_mov_11';
% frcPath = 'H:\tlukes\data_Hendrik_livecells\for_tomas1\fourier ring correlation\20150722_mov_3';

% stack = load_tifFile([imagePath, filesep,imageName,'.tif'],[],[]);
% load([imagePath, filesep,imageName,'.mat']);

dip_initialise
tic
% settings.orders = 1:5;
% settings.bcgsub = 0;
% sample = 'mEos2';
% sample = 'psCFP2';
% note = 'bcgsub1x_adaplin';
% note = '_adaplin';
%%
cut = 1;

for mm = settings.orders%3:4
    order = mm;
    if order == 1
        stack = sofi_c{order};
    else
        stack = sofi{order}(cut:end,cut:end,:);
    end
    [sy,sx,frames] = size(stack); 

%     if strcmp(sample, 'psCFP2')
%     % x74 y10 w20 h10
%     stack(5*mm+1:5*mm+5*mm,37*mm+1:37*mm+10*mm,:) = zeros(5*mm,10*mm,frames);
%     end
    
    %%% Calculate FRC
        % First do not forget to initialize dipimage package by writing 
        % "dip_initialise" into command window

        fprintf('\n -- computing FRC curve--\n')
        nss = size(stack,3);
        sr = 1;%mm-2;
        io = sr*order;
%         pixelsize = 104.8; %nm
        superzoom = io;
        roix = 1:sx;
        roiy = 1:sy;
        numSubSets = 1:nss;
        numTests = 10;
        thresh = 0.0;
        
%         settings.roix = 40*io:40*io + 100*io; 
%         settings.roiy = settings.roix;
%         stack = stack(roiy,roix,:);
    count = 1; 
    w = 4;
    step = 2;
    for kk = nss%4:2:nss %numSubSets
                disp(['load block number',num2str(kk)])
                in1 = mean(stack(:,:,1:2:kk),3); 
                in2 = mean(stack(:,:,2:2:kk),3);
%                 in1 = mean(stack(:,:,kk-4+1:2:kk),3); 
%                 in2 = mean(stack(:,:,kk+2-4:2:kk),3);

%                 in1 = mean(stack(:,:,kk-step+1:2:kk+step),3); 
%                 in2 = mean(stack(:,:,kk-step+2:2:kk+step),3);
%                 figure, 
%                 subplot(121);
%                 imshow(in1,[]);
%                 subplot(122);
%                 imshow(in2,[]);
%                 in1in2 = mean(cat(3,in1,in2),3);

                sc_in1 = max(in1(:));
                sc_in2 = max(in2(:));
                in1 = in1./sc_in1;
                in2 = in2./sc_in2;
%                 in1 = imadjust(in1,[thresh 1],[0 1]);
%                 in2 = imadjust(in2,[thresh 1],[0 1]);
                if settings.bcgsub > 0 && order > 1
                bcg1 = median(in1(:));
                bcg2 = median(in2(:));
                in1 = imadjust(in1,[settings.bcgsub*bcg1 1],[0 1]);
                in2 = imadjust(in2,[settings.bcgsub*bcg2 1],[0 1]);
                end
%                 in1 = imresize(in1,sr,'bilinear');
%                 in2 = imresize(in2,sr,'bilinear');

%                 in1m = in1(1:sy*sr,1:sx*sr,:);
%                 in2m = in2(1:sy*sr,1:sx*sr,:);

%                 figure, 
%                 subplot(121);
%                 imshow(in1,[]);
%                 subplot(122);
%                 imshow(in2,[]);

            
                sector = [0,1*pi/6;
                          1*pi/6, 2*pi/6;
                          2*pi/6, 3*pi/6;
                          3*pi/6, 4*pi/6;
                          4*pi/6, 5*pi/6;
                          5*pi/6, 6*pi/6];
                      
%                       sector = [0,1*pi/6;
%                           1*pi/6, 2*pi/6;
%                           2*pi/6, 3*pi/6;
%                           3*pi/6, 4*pi/6;
%                           4*pi/6, 5*pi/6;
%                           5*pi/6, 6*pi/6;
%                           6*pi/6, 7*pi/6;
%                           7*pi/6, 8*pi/6;
%                           8*pi/6, 9*pi/6;
%                           9*pi/6, 10*pi/6;
%                           10*pi/6, 11*pi/6;
%                           11*pi/6, 12*pi/6];
%                 sector = [-2*pi/6,0,0,2*pi/6];    

                for isector = 1:size(sector,1)+1     
                if isector==size(sector,1)+1 
                    [resolution, resolution_h, resolution_l,frc_curve] = frcsec(in1,in2,[]);
                else
                    [resolution, resolution_h, resolution_l,frc_curve] = frcsec(in1,in2,sector(isector,:));  
                end
%                 in1 = dip_image(in1m);
%                 in2 = dip_image(in2m);
% 
%                 sz = imsize(in1);
% 
%                 % Calculate FRC curve
% 
%                 % Compute mask in x-direction
%                 nfac = 8;                                                   % Image width / Width of edge region
%                 x_im = xx(sz(1),sz(2))/sz(1);
%                 mask = 0.5-0.5*cos(pi*nfac*x_im);          
%                 mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;
% 
%                 % Check that input images are square and mask
%                 if sz(1) == sz(2)
%                     mask = mask*rot90(mask);
% 
%                     % Mask input images
%                     in1 = mask*in1;
%                     in2 = mask*in2;
%                 else
%                     warning('frc:nonsquare','Images are not square.');
% 
%                     % Compute mask in y-direction
%                     y_im = yy(sz(1),sz(2))/sz(2);
%                     mask_y = 0.5-0.5*cos(pi*nfac*y_im);          
%                     mask_y(abs(y_im)<((nfac-2)/(nfac*2))) = 1;
%                     mask = mask*mask_y;
%                     clear mask_y
% 
%                     % Mask input images
%                     in1 = mask*in1;
%                     in2 = mask*in2;
% 
%                     % Make images square through zero padding
%                     in1 = extend(in1,[max(sz) max(sz)]);
%                     in2 = extend(in2,[max(sz) max(sz)]);
%                 end
% 
%                 % Fourier transform input images
%                 in1 = ft(in1);
%                 in2 = ft(in2);
% 
%                 % Compute fourier ring correlation curve
%                 frc_num = real(radialsum(in1.*conj(in2)));                                  % Numerator
%                 in1 = abs(in1).^2;
%                 in2 = abs(in2).^2;
%                 frc_denom = sqrt(abs(radialsum(in1).*radialsum(in2)));                      % Denominator
%                 frc_out = double(frc_num)./double(frc_denom);                               % FRC
%                 frc_out(isnan(frc_out)) = 0;                                                % Remove NaNs
% 
%                 %
%                 % Calculate the resolution
%                 sz = max(imsize(in1));
% 
%                 % resolution = frctoresolution(frc_out,sz);
% %                 frc_out(end-10:end) = 0;
%                 [resolution, resolution_h, resolution_l] = frctoresolution(frc_out,sz);
                %
                resMid(isector) = resolution*(pixelsize./superzoom);
                resHigh(isector) = resolution_h*(pixelsize./superzoom);
                resLow(isector)  = resolution_l*(pixelsize./superzoom);
                frc_curves(isector,:) = frc_curve;
%                 resolutionAll = mean(resolution_mid)
                end
    % [~,frc_curve] = postoresolution(coords, szx, superzoom); 

    % in1in2 = mean(cat(3,in1,in2),3);
%     in1in2 = imadjust(in1in2./max(in1in2(:)),[thresh 1],[0 1]);

    %% calc mean frc and show results for kk-th subset

%         figure('Visible','On'), 
%     %     in1in2 = deconvlucy(in1in2,fspecial('gaussian',[29 29],1.3),10);
%     %      in1in2 =  in1in2.^(1/4);
%         subplot(221);imshow(in1in2,[]);title('analysed region');hold on; 
% 
%             hold on;
%             y_vals = roiy;
%             x_vals = roix;
% 
%             plot(x_vals,ones(1,length(x_vals))*min(y_vals),'-r');
%             plot(x_vals,ones(1,length(x_vals))*max(y_vals),'-r');
%             plot(ones(1,length(y_vals))*min(x_vals),y_vals,'-r');
%             plot(ones(1,length(y_vals))*max(x_vals),y_vals,'-r');
%             hold off;
% 
%         subplot(222);imshow(in2m,[]);title('analysed region in detail');
% 
%         subplot(223)
%             % Show frc curve
%             frc_curve = frc_out;
%             qmax = 0.5/(pixelsize/superzoom);
%             plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
%             xlim([0,qmax])
%             hold on
%             plot([0 qmax],[1/7 1/7],'r-');
%             plot([0 qmax],[0 0],'k--'); hold off
%             xlabel('spatial frequency (nm^{-1})')
%             ylabel('FRC')
%             title('Fourier Ring Correlation curve')
% 
%         % Show measured resolution
%         subplot(224)
%             res = [resolution_high;resolution_mid;resolution_low];
        %     res(:,res(1,:)==0) = [];
%             res_mean = median(res,2);
% %             res_std = std(res,[],2);
%             res1 = res(1,:);
%             res2 = res(2,:);
%             res3 = res(3,:);
%             res1(isnan(res1)) = [];
%             res2(isnan(res2)) = [];
%             res3(isnan(res3)) = [];
            
        %     res1(res(1,:)>1.5*res_mean(1)) = [];
%             res2(res(2,:)>1.5*res_mean(2)) = [];
%             res3(res(3,:)>1.5*res_mean(3)) = [];
% 
%             res1(res1==0) = [];

%             resHigh = mean(res1);
%             resMid = mean(res2);
%             resLow = mean(res3);
% 
%             h1 = plot(res(1,:),'+-.g');hold on;
%             h2 = plot(res(2,:),'+--b');
%             h3 = plot(res(3,:),'+--m');
%             ylabel('resolution in [nm]');
%             xlabel('measurement number [-]');
%             legend('Location','SouthOutside','resolution high', 'resolution mid', 'resolution low');
%             errorbar(1:length(res), res(1,:), ones(1,length(res))*res_std(1)/2, ones(1,length(res))*res_std(1)/2, 'r', 'Marker', 'none', 'LineStyle', 'none' );
%             errorbar(1:length(res), res(2,:), ones(1,length(res))*res_std(2)/2, ones(1,length(res))*res_std(2)/2, 'r', 'Marker', 'none', 'LineStyle', 'none' );
%             errorbar(1:length(res), res(3,:), ones(1,length(res))*res_std(3)/2, ones(1,length(res))*res_std(3)/2, 'r', 'Marker', 'none', 'LineStyle', 'none' );
%             title(['resolution high, mid, low :', num2str(resHigh),'  ', num2str(resMid),'  ',num2str(resLow)]); 

            % save the current figure
%             fh = gcf;saveFigure(fh,outputPath,'FRC_resolution','png',1); 

            results(count,mm).resHigh = resHigh;
            results(count,mm).resMid = resMid;
            results(count,mm).resLow = resLow;
            results(count,mm).frc_curves = frc_curves;
            count = count+1;
            clear frc_curves;
%             results(kk,mm).resHigh_std = std(res1);
%             results(kk,mm).resMid_std = std(res2);
%             results(kk,mm).resLow_std = std(res3);
            
    end
end
toc


%%
%         figure, 
%         linecolor= {'b.:','g.:','r.:','m.:','k.:'};  
%         count = 1;
%         for mm = settings.orders
% %             res_std = [results(numSubSets,mm).resMid_std];
%             res_mean = [results(numSubSets,mm).resMid];
%             
%     %         h(mm) = plot([results(numSubSets,mm).resMid],[linecolor{mm-2},'.:']);
%             hold all;
%             
%     %         set(gca,'XTickLabel',['0';' ';'1';' ';'2';' ';'3';' ';'4'])
%             xframes = 0:2000:20000;
%     %         xlim([1 20]);
%             plot(numSubSets, res_mean, linecolor{count});
%             legendInfo{count} = ['order ',num2str(mm)];
% %             res_std = [results(numSubSets,mm).resLow_std];
% %             res_mean = [results(numSubSets,mm).resHigh];
%             
% %             errorbar(numSubSets, res_mean, res_std./2, res_std./2, linecolor{count});
%             
% %             res_mean = [results(numSubSets,mm).resLow];
%             
% %             errorbar(numSubSets, res_mean, res_std./2, res_std./2, linecolor{count});
%             count = count+1;
%         end
%         hold off
%         grid on;
%         ylabel('FRC resolution [nm]');
%         xlabel('Number of frames [-]');
%         set(gca,'XLim',[0 20])
%         set(gca,'YLim',[0 800]);
%         
%         set(gca,'XTick',[0:2:20])
%         set(gca,'YTick',[0:100:800])
%         set(gca,'XTickLabel',xframes);
%         
%         %         title('psCFP2');
%         title(sample);
%         
% %         legend([h{1},h{2},h{3}],legendInfo{1},legendInfo{2},legendInfo{3});
% %         hleg = legend(legendInfo);
%         figname = ['_FRC_bcgsub_',num2str(settings.bcgsub),'_',note];
%         fh = gcf;saveFigure(fh,outputPath,figname,'fig'); 
%         fh = gcf;saveFigure(fh,outputPath,figname,'png');
%         
% %% 
% 
%         linecolor= {'b.:','g.:','r.:','m.:','k.:'};  
%         count = 1;
%         for mm = settings.orders
%           figure, 
% %             res_std = [results(numSubSets,mm).resMid_std];
%             res_mean = [results(numSubSets,mm).resMid];
%             
%     %         h(mm) = plot([results(numSubSets,mm).resMid],[linecolor{mm-2},'.:']);
%             hold all;
%             
%     %         set(gca,'XTickLabel',['0';' ';'1';' ';'2';' ';'3';' ';'4'])
%             xframes = 0:2000:20000;
%     %         xlim([1 20]);
%             plot(numSubSets, res_mean,linecolor{count});
%             legendInfo{count} = ['Mid estimate, order ',num2str(mm)];
% %             res_std = [results(numSubSets,mm).resLow_std];
%             res_mean = [results(numSubSets,mm).resLow];
%             
%             plot(numSubSets, res_mean,linecolor{count+1});
%             legendInfo{count+1} = ['High estimate,order ',num2str(mm)];
% %             res_std = [results(numSubSets,mm).resHigh_std];
%             res_mean = [results(numSubSets,mm).resHigh];
%             
%             plot(numSubSets, res_mean,linecolor{count+2});
%             legendInfo{count+2} = ['Low estimate,order ',num2str(mm)];
% %             count = count+1;
%             hold off
%             ylabel('FRC resolution [nm]');
%             xlabel('Number of frames [-]');
%             set(gca,'XLim',[0 20]);
%             set(gca,'YLim',[0 800]);
%             set(gca,'YTick',[0:100:800]);
%             set(gca,'XTick',[0:2:20]);
%             set(gca,'XTickLabel',xframes);
%             title(sample);
%             legend(legendInfo);
%             grid on;
%             % save the current figure
%             
%             figname = ['_FRC_io',num2str(mm),'_bcgsub_',num2str(settings.bcgsub),'_',note];
%             fh = gcf;saveFigure(fh,outputPath,figname,'fig'); 
%             fh = gcf;saveFigure(fh,outputPath,figname,'png'); 
%         end
       
    
%%
save([outputPath,filesep,'FRC_bcgsub_',num2str(settings.bcgsub),'_',settings.note,'.mat'],'results','settings','-v7.3');
   
%% plot FRC
clear frc;
clear legendinfo;

io=2;
for ii = 1:size(results,1)
    for jj=4
   frc(ii,1)=results(ii,io).resMid(jj);
   legendinfo{1} = ['Sector: ',num2str(jj)];
    end
    
end

figure, 
plot(1:length(frc),frc,'x--');
xlabel('Frame of SR sequence');
ylabel('sFRC');
legend(legendinfo);

ylim([100 300]);
% xlim([0.5 4.5]);
set(gca,'XTick',[1:7])
% xframes = [1, 2, 3, 4,];
% set(gca,'XTickLabel',xframes);
text(1:length(frc),frc,num2str(round(10*frc)./10),'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
    
       
