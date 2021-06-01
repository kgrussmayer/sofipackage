% 2D SOFI - batch process all input files
% Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2012a

clear all;
% close all;clc;

c1=now;
warning('off','all')

% names of all the files to be processed
% % fnames = {'20mW_DU897_BV_2342',...
% %           '20mW_DU897_BV_2343',...
% %           '20mW_DU897_BV_2345',...
% %           '20mWEMG50_DU897_BV_2347',...
% %           '20mWEMG50_DU897_BV_2348',...
% %           '20mWEMG50_DU897_BV_2349',...
% %           '100mW_DU897_BV_2341'};
% % fnames = {'130806_PALM_LATwtFL-mEos2_EEZetaFL-PSCFP2_dual--COS7_fixed_PSCFP2_ser04_crop',...
% % '130813_PALM_CD4wt-mEos2_LATwtFL-PSCFP2_dual--COS7_fixed_mEos2_ser03_crop',...
% % '130813_PALM_CD4wt-mEos2_LATwtFL-PSCFP2_dual--COS7_fixed_PSCFP2_ser03_crop',...
% % 'PALM_130410_CD4wt-mEos2_fixed_Jurkat_ser08',...
% % 'PALM_130425_LATFL-mEos2_fixed_HeLa_ser06_crop',...
% % 'PALM_130425_LATFL-mEos2_fixed_HeLa_ser07_crop',...
% % 'PALM_130502_LAT-FL-mEos2_fixed_Jurkat_rest_ser03_crop',...
% % 'PALM_130515_EEzetaFL-PSCFP2_HeLa_fixed_ser01_crop',...
% % 'PALM_140320_pXJCD4dCTPSCFP2_Jur_fix_ser07',...
% % 'PALM_CD4wt-PAmCh_COS7_fixed_ser7'};
% % fnames = {'MitoOrange_40nM_532_ND1_405_2mW_10ms_DU897_BV_2393',...
% %     'MitoOrange_40nM_532_ND1_405_2mW_10ms_DU897_BV_2394',...
% %     'MitoOrange_40nM_532_ND1_405_5mW_10ms_DU897_BV_2391',...
% %     'MitoOrange_40nM_532_ND1_405_5mW_10ms_DU897_BV_2392',...
% %     'MitoOrange_40nM_532_ND1_405_10mW_10ms_DU897_BV_2389',...
% %     'MitoOrange_40nM_532_ND1_405_10mW_10ms_DU897_BV_2390'};
% % fnames = {'MEA20mMOSS30ms3000imgsND06800mWGreen_DU897_BV_0917',...
% %     'MEA50mMOSS30ms3000imgs800mWND06Green_DU897_BV_0919'};
% % fnames = {'3T3_Lyn_mOrange_30ms_532nm_ND2.5_DU897_BV_2368',...
% %           '3T3_Lyn_mOrange_30ms_532nm_ND2.5_DU897_BV_2369',...
% %           '3T3_Lyn_mOrange_30ms_532nm_ND3_DU897_BV_2366',...
% %           '3T3_Lyn_mOrange_30ms_532nm_ND3_DU897_BV_2367',...
% %           '3T3_Lyn_mOrange_30ms_532nm_ND3_DU897_BV_2367_X2'};
% % fnames = {'30msMEA20mMOSSGreenND03UV20mW_DU897_BV_1790',...
% %           '30msMEA20mMOSSGreenND03UV20mW_DU897_BV_1791',...
% %           '30msMEA20mMOSSGreenND03UV20mW_DU897_BV_1792',...
% %           '30msMEA20mMOSSGreenND03UV20mW_DU897_BV_1793',...
% %           '30msMEA20mMOSSGreenND03UV20mW_DU897_BV_1794',...
% %           '30msMEA20mMOSSRedUV20mW_DU897_BV_1789'};   
% % fnames = {'30msMEA50mMOSSGreen_ND1_DU897_BV_1780'}; 
% % fnames = {'30msMEA50mMOSSRedUV5mW_DU897_BV_1772',...
% %           '30msMEA50mMOSSRedUV5mW_DU897_BV_1773',...
% %           '30msMEA50mMOSSRedUV1mW_DU897_BV_1774',...
% %           '30msMEA50mMOSSRedUV1mW_DU897_BV_1775',...
% %           '30msMEA50mMOSSRedUV1mW_DU897_BV_1776'};     
% % fnames = {'Hela_Pmi_5uM_HBSS_635nm_FP_15ms_DU897_BV_2397',...
% %           'Hela_Pmi_5uM_HBSS_635nm_FP_ND0_15ms_DU897_BV_2399',...
% %           'Hela_Pmi_5uM_HBSS_635nm_FP_ND04_15ms_DU897_BV_2398'};
% % fnames = {'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0_15ms_DU897_BV_2400';
% %           'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0_15ms_DU897_BV_2401'};   
% % fnames = {'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0_405nm_5mW_15ms_DU897_BV_2402',...
% % 'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0_405nm_1mW_15ms_DU897_BV_2403',...
% % 'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0.4_405nm_1mW_15ms_DU897_BV_2404'};    
% % fnames = {'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0_405nm_1mW_15ms_DU897_BV_2405',...
% % 'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2406',...
% % 'Hela_Pmi_1.2uM_Leib_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2407'};     
% % fnames = {'Hela_Pmi_1.2uM_L0ib_buf1_635nm_FP_ND0.3_405nm_1mW_7.5ms_DU897_BV_2412',...
% %   'Hela_Pmi_1.2uM_L0ib_buf1_635nm_FP_ND0.3_405nm_1mW_7.5ms_DU897_BV_2413',...
% %   'Hela_Pmi_1.2uM_L0ib_buf1_635nm_FP_ND0_405nm_1mW_7.5ms_DU897_BV_2410',...
% %   'Hela_Pmi_1.2uM_L0ib_buf1_635nm_FP_ND0_405nm_1mW_7.5ms_DU897_BV_2411',...
% %   'Hela_Pmi_1.2uM_Leib_buf1_635nm_FP_ND0_405nm_1mW_15ms_DU897_BV_2409'};
% % fnames = {'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2414',...
% %        'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2415'};                 
% % fnames = {         
% %    'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0_405nm_1mW_7.5ms_DU897_BV_2416',...
% %    'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0_405nm_1mW_7.5ms_DU897_BV_2417',...
% %    'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0_405nm_1mW_7.5ms_DU897_BV_2418',...
% %    'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0_405nm_1mW_7.5ms_DU897_BV_2419'};
% % fnames = {
% %    'Hela_Pmi_1.2uM_Leib_buf2_635nm_FP_ND0_405nm_3mW_7.5ms_DU897_BV_2420'};
% % fnames = {'Hela_Pmi_1.2uM_Leib_buf3_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2420',...
% %    'Hela_Pmi_1.2uM_Leib_buf3_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2421',...
% %    'Hela_Pmi_1.2uM_Leib_buf3_635nm_FP_ND0_405nm_1mW_15ms_DU897_BV_2422'};
% % fnames = {'Hela_Pmi_1.2uM_Leib_buf4_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2423',...
% %    'Hela_Pmi_1.2uM_Leib_buf4_635nm_FP_ND0.3_405nm_1mW_15ms_DU897_BV_2424',...
% %    'Hela_Pmi_1.2uM_Leib_buf4_635nm_FP_ND0_405nm_1mW_15ms_DU897_BV_2425'};
% % fnames = {'GreenND02_405ND01_365ND11_DU897_BV_1544','GreenND02_405ND05_365ND09_DU897_BV_1561','GreenND02_405ND13_365ND07_DU897_BV_1592',...
% %    'GreenND02_405ND17_365ND05_Luc247_MONO_0468'};
% % fnames = {'GreenND02_405ND08_365ND08_DU897_BV_1580','GreenND02_405ND11_365ND07_DU897_BV_1584'};
% % fnames = {'MEF_Lifeact_meos2_30ms_488nm_20mW_DU897_BV_2387','MEF_Lifeact_meos2_30ms_488nm_5mW_DU897_BV_2386',...
% %    'MEF_Lifeact_meos2_30ms_488nm_5mW_DU897_BV_2385','MEF_Lifeact_meos2_30ms_488nm_10mW_DU897_BV_2384','MEF_Lifeact_meos2_30ms_488nm_10mW_DU897_BV_2383',...
% %    'MEF_Lifeact_meos2_30ms_488nm_25mW_DU897_BV_2382','MEF_Lifeact_meos2_30ms_488nm_25mW_DU897_BV_2381',...
% %    'MEF_Lifeact_meos2_30ms_488nm_25mW_DU897_BV_2380','MEF_Lifeact_meos2_30ms_488nm_20mW_DU897_BV_2379'};
% % fnames = {'GreenND02_405ND01_365ND11_DU897_BV_1544','GreenND02_405ND05_365ND09_DU897_BV_1561','GreenND02_405ND11_365ND07_DU897_BV_1584',...
% %    'GreenND02_405ND17_365ND05_DU897_BV_1566'};   

% fnames = {'mol1','mol2','mol3','mol4'};
% fnames = {'mol1','mol2','mol4'};
% fnames = {'20150619_fixed_cell_meos2_paxillin_fiducial_mov_10','20150619_fixed_cell_meos2_paxillin_fiducial_mov_11',...
%     '20150619_living_cell_meos2_paxillin_fiducial_mov_3',...
%     '20150722_mov_1','20150722_mov_2','20150722_mov_3'};

% fnames = {'20150619_fixed_cell_meos2_paxillin_fiducial_mov_10','20150619_fixed_cell_meos2_paxillin_fiducial_mov_11'};
% fnames = {'20150619_fixed_cell_meos2_paxillin_fiducial_mov_11'};

% fnames = {'20150722_mov_1','20150722_mov_3','20150722_mov_2'};
% fnames = {'20150722_mov_2b','20150722_mov_3b'};
% fnames = {'20150722_mov_2b'};
fnames = {'20150722_mov_3b'};
% fnames = {'20150902_mov_2'};
% fnames = {'20150902_mov_5'};
%     '20150619_living_cell_meos2_paxillin_no_fiducial_mov_1',...
%     '20150619_living_cell_meos2_paxillin_no_fiducial_mov_2'};

% folder that contains files to be processed
% % imagePath='H:\tlukes\data_STORM';
% % imagePath='H:\Sharipov\Measurements\20150610_Fps_Arik';
% % imagePath='H:\tlukes\data_HOF_JHI_sample';
% % imagePath='H:\Sharipov\Measurements\20150626_Hela_Arik\Mitotracker';
% % imagePath='I:\sharipov\20130814_HepaBlinking\Atto532_postfixation';
% % imagePath='I:\sharipov\20140521_Hela_TOM20(Tris_buffer)';
% % 
% % 
% % imagePath='I:\sharipov\20140519_Hela_TOM20'; %%TODO
% % imagePath='I:\sharipov\20140519_Hela_DM1A';
% % imagePath='H:\Sharipov\Measurements\20150630_Hela_mitostaining_Sandra';
% % imagePath='H:\Sharipov\Measurements\20150630_Hela_mitostaining_Sandra\Buffer4';
% % imagePath='H:\Sharipov\Measurements\!Imgs_onratios_for_Tomas';
% % imagePath='H:\Sharipov\Measurements\20150617_Fps_Arik';
% % imagePath='H:\Sharipov\Measurements\20150611_Hela_3T3_Lyn_mOrange\3T3_Lyn_mOrange';

% imagePath='H:\tlukes\data_DNA sequencing\T7';


fnames = {'MEF_10ms532nm_FP_ND2_405_50mw_ND1_DU897_BV_2576'};


imagePath='H:\tlukes\data_Hendrik_livecells';
imagePath='H:\tlukes\data_Hendrik_livecells\movies';


% I/O settings 
settings.io.outputpath = [imagePath,filesep,'SNR_seq500_fwhm3_mov3_roi_2'];% output folder for results
% settings.io.outputpath = [imagePath,filesep,'results_subseq500_fwhm3_ksofi_4_ksofialg'];% output folder for results
% settings.io.outputpath = ['H:\tlukes\results2D\onratioTest',filesep,'results_subseq500_fwhm3'];
settings.io.bits = 16; % number of bits of the output tif file {8,16}
% settings.io.roi = [201:300;201:300]; % region of interest to process - keep empty [] if the whole image should be used
settings.io.roi = []; % process the whole image
% settings.io.roisx = { 1+40:256-36,1+8:256-12,1+11:256-14,...
%                       1+10:257-13,1+12:257-54,1+11:257-13};
% 
% settings.io.roisy = {1+41:301-7,1+150:301,1+60:301-72,...
%                     1:301-30,1:301-14,1:301};
%                 
% settings.io.roisx = { 1+8:256-12,1+11:256-14,...
%                       1+10:257-13,1+12:257-54,1+11:257-13};
% 
% settings.io.roisy = {1+150:301,1+60:301-72,...
%                     1:301-30,1:301-14,1:301};
                

% settings.io.roisx = {94:244}; % mov11 roi 1
% settings.io.roisy = {151:301};

% settings.io.roisx = {170:210}; % mov11 roi 2
% settings.io.roisy = {201:241};
% 
% settings.io.roisx = {100:140}; % mov11 roi 3
% settings.io.roisy = {196:236};
% 
% settings.io.roisx = {115:155}; % mov11 roi 4
% settings.io.roisy = {253:293};

% settings.io.roisx = {32:182}; % mov3b roi 1
% settings.io.roisy = {141:291};
% 
% settings.io.roisx = {40:80}; % mov3b roi 2
% settings.io.roisy = {232:272};

settings.io.roisx = {353:468}; % live movie 2576 roi 3
settings.io.roisy = {6:121};

% settings.io.roisx = {83:123}; % mov3b roi 3
% settings.io.roisy = {210:250};

% settings.io.roisx = {141:181}; % mov3b roi 4
% settings.io.roisy = {187:227};

% settings.io.roisx = {44:230}; % 20150902_mov2 roi 1
% settings.io.roisy = {71:257};

% settings.io.roisx = {34:220}; % 20150902_mov_5 roi 1
% settings.io.roisy = {95:281};

% settings.io.roisx = {1+12:257-54,1+11:257-13};
% settings.io.roisy = {1:301-14,1:301};
% 
% settings.io.roisy = {1:301};
% settings.io.roisx = {1+11:257-13};


settings.io.blcor = 0; % bleaching correction off/on {0,1}
settings.io.dcor = 1;% turn on or of drift correction, if on - specify path to calib file


% Cumulant calculation settings
settings.sys.orders = 1:2;
settings.sys.wsize = 500;
settings.sys.sub = []; % evaluate only first n frames (for quick preview)
settings.sys.start = 1;
settings.sys.jk = 1; % turn on/off the Jacknife SNR estimation
settings.sys.block = 1;

% Deconvolution settings
settings.dec.fwhm = 4;
settings.dec.iter = 10;
settings.dec.gammas = [1, 0.5, 0.55, 0.4];
% Molecular parameters
settings.molpar.thresh = 0.2;
settings.molpar.run = 0;

%%
for ii = 1:numel(fnames)
    c(ii)=now;
    disp(['Processing file number: ',num2str(ii),' from ',num2str(numel(fnames))]);
    settings.io.imageName = fnames{ii}; 
    settings.io.imageFile = [imagePath, filesep,settings.io.imageName];
    settings.io.roit = {settings.io.roisx{ii},settings.io.roisy{ii}}; 
%     [sofi,sofi_c,settings] = SOFI2D_process(settings);
    [sofi,sofi_c,settings,snrs,vars,bias] = SOFI2D_process_jk(settings); 
%     [sofi,sofi_c,settings,snrs,vars,bias,snrs2,vars2] = SOFI2D_process_jkseq_02(settings);
%     [sofi,sofi_c,settings,snrs,vars,bias] = SOFI2D_process_jkseq_05(settings); 
    
%     save([settings.io.outputpath,filesep,settings.io.imageName,'.mat'],'sofi','sofi_c')
%     save([settings.io.outputpath,filesep,settings.io.imageName,'_snr_roi2_500_b1_all_o1-3_dcon_bcoff_m3.mat'],'snrs','vars','bias')
end

c2=now;

disp(['Total time [hours]:',num2str((c2-c1)*24)])
%% show results snr
% count = 1;
% figure,
% for ii = 3:max(settings.sys.orders)
% %     plot(10*log10(squeeze(mean(mean(sofi_snr{ii},1),2)))); hold on;
%     plot(squeeze(mean(mean(sofi_snr{ii},1),2))); hold on;
%     legendInfo{count} = ['SOFI order ',num2str(ii)];
%     count = count +1;
% end
% legend(legendInfo)
% title('Mean SNR')
% % ylabel('Mean SNR [dB]');
% xlabel('Number of a subsequence [-]')

%% show results
count = 1;
figure('Position',[200,400,1300,500]),

xaxis = (1:size(snrs,2)).*settings.sys.wsize;
for io = settings.sys.orders
    for pp = 1:size(snrs,3)
        subplot(1,2,io)
        for ii = 1:size(snrs,2)
            temp =snrs{io,ii,pp};
            temp(isnan(temp)) = 0;
%         plot_snr(ii) = mean2(temp(10:end-10,10:end-10));
        plot_snr(ii) = mean2(temp);
%         plot_snr(ii) = mean2(vars{io,ii,pp});
%         plot_snr2(ii) = mean2(snrs2{io,ii,pp});
        end
        plot_snr(isinf(plot_snr)) = 0;
%         if io ==1
        plot(xaxis,plot_snr,'Marker','x','LineStyle','--');hold on;
%         else
%         plot(xaxis,100*1./plot_snr,'b--+');hold on;
%         end
%         pause
        legendInfo{count} = ['snr order: ',num2str(io),' block size ',num2str(settings.sys.block.*104.8),'nm'];
%         plot(xaxis,plot_snr2,'r--');hold on;
%         legendInfo{count+1} = ['method2 snr order: ',num2str(io),' block size ',num2str(0)];
        xlabel('Number of frames');
        ylabel('SNR [-]');
        title(legendInfo{count});
        count = count+1;
    end
end
% xlabel('Frames');
% ylabel('SNR [-]');
% legend(legendInfo);

%% 
figure, 
subplot(131);
imshow(sqrt(vars{2,1,pp}),[]);colormap('jet');colorbar;
subplot(132);
temp = binStack(abs(sofi_c{2}(:,:,1)),2).^settings.dec.gammas(2);
imshow(temp,[]);colormap('jet');colorbar;
subplot(133);
imshow(temp./sqrt(vars{io,1,pp}),[]);colormap('jet');colorbar;

mean2(temp./sqrt(vars{io,1,pp}))
%% Save settings and processing info
settings.io.proctime = (c2-c1)*24;

save([settings.io.outputpath,filesep,'config.mat'],'settings'); % write to a mat file
descfile_sofi3D(settings); % write to a text file

%% JK1

figure, 
subplot(131);
imshow(snrs{2}(:,:,1),[]);colormap('jet');colorbar;
subplot(132);
imshow(vars{2}(:,:,1),[]);colormap('jet');colorbar;
subplot(133);
imshow(bias{2}(:,:,1),[]);colormap('jet');colorbar;



