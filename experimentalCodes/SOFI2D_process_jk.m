% 2D SOFI processing 
% Calculate SOFI cumulants, perform flattening, deconvolution and
% linearization accoring to bSOFI approach
% Tomas Lukes, tomas.lukes@epfl.ch
% created: 15.4.2015

function [sofi,sofi_c,settings,sofi_snr,sofi_vars,sofi_bias] = SOFI2D_process_jk(settings) 


    %% Load image stack

    stack = load_tifFile([settings.io.imageFile,'.tif'],[],settings.io.roi);
%     [stack,frames] = load_bigtif([settings.io.imageFile,'.tif']);
    [sy,sx,frames] = size(stack); 
%     frames = 1000;
    %% Bleaching correction
    % simple gaussian fitting to the curve given by the decay of average intensity for each frame 
    if settings.io.blcor == 1
        mitrace=squeeze(mean(mean(stack,1),2));
        mitrace = mitrace./max(mitrace);
        s = fitoptions('Method','NonlinearLeastSquares',...
                       'Lower',[0,0,0],...
                       'Upper',[Inf,Inf,Inf],...
                       'Startpoint',[1 1 1]);
        f = fittype('a*exp(-x/b)+c','options',s);

        [c2,gof2] = fit((0:frames-1)',mitrace,f);
        a = c2.a;
        b = c2.b;
        c = c2.c;

        figure;x=(0:frames-1);plot(x',mitrace);hold on;plot(x,a*exp(-x/b)+c,'k');

        for ii = 1:frames
            stack(:,:,ii) = stack(:,:,ii)/(a*exp(-(ii-1)/b)+c);
        end
    end
    
    %% Drift correction 
    % drift correction based on tracking fiducial markers
    disp(frames);
    if settings.io.dcor == 1
    load([settings.io.imageFile,'_drift_corr.mat'])
%     drift = drift_corr';
    drift = drift_corr;
    %%% correct drift
    pxsize = 104.8; % relative pixel size in nm
    [sy,sx,frames] = size(stack);
    [xx, yy]=meshgrid(1:sx,1:sy);

    for k=1:size(drift,1)
        disp(['dcor: ',num2str(k)]);
        img = stack(:,:,k+1);

        %%% total drift
%         drift_subseq = sum(drift(1:k,:),1);

        shiftxy = drift(k,2:3);
        shiftxy = shiftxy./(pxsize);

        imReg=interp2(xx,yy,double(img),xx+shiftxy(1),yy+shiftxy(2),'linear');
        imReg(isnan(imReg)) = 0;
        stack(:,:,k+1)=uint16(imReg);
    end
    
    if isfield(settings.io,'roit')
        roi = settings.io.roit;
        stack = stack(roi{2},roi{1},:);
    end
    
    end
    
%     [sy,sx,frames] = size(stack);


    %% take only substack
%     stack = stack(:,:,[6500:7000,9500:10000]); % mov11 sub
%     stack = stack(:,:,[4000:4500,13500:14000]); % mov3b sub
  
%     [sy,sx,frames] = size(stack);
    %% Calculate cumulants, perform flattening

    count =1;
    for jj = settings.sys.start:settings.sys.wsize:floor((frames-settings.sys.start+1)/settings.sys.wsize)*settings.sys.wsize
        substack = stack(:,:,jj:jj+settings.sys.wsize-1);

        tic
        if settings.sys.jk == 1
            [c,grid,bias,snrs,vars]=sofiCumulantsSNR(substack,1,[],[],settings.sys.orders);
        else 
            [c,grid]=sofiCumulants2D(substack,1,[],[],settings.sys.orders); % cross cumulants, zero time lag 
    %     [cf,grid]=fourierCumulants(substack,1,[],[],orders);
%         ck{2} = specSofi2D(substack);
%         ck{2} = ksofi(substack);
        end
        toc

        c=sofiAllFlatten(c,settings.sys.orders);

        for io=settings.sys.orders
            sofi_c{io}(:,:,count)=c{io};
            if settings.sys.jk ==1;
            sofi_snr{io}(:,:,count)=snrs{io};
            sofi_vars{io}(:,:,count)=vars{io};
            sofi_bias{io}(:,:,count)=bias{io};
            end
%             sofi_ck{io}(:,:,count)=ck{io};
    %         sofi_cf{io}(:,:,count)=cf{io};
        end
        count = count+1;
    end

    %% Deconvolution and postprocessing
%     num = 500;
%     [gcor,rho] = kest_v2(sofi_c,num);
%     gammas = gcor./repmat((1:settings.sys.orders(end))',1,size(gcor,2));
%     gammas(1,:)=1;

%     [sofi_lin,imgpar2,imgpar3,gamcor]=sofiLinearize3(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,gammas);

    %     sofi=sofiLinearize_my(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter);
    sofi=sofiLinearize(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter);
%     sofif=sofiLinearize(sofi_ck,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter);
    
%     sofi2=sofiLinearize_my(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter);
%     iter = 100;
%     gammas = [1 0.5 0.45 0.35];
%     lincoeffs = [1 0.5 0.35 0.35];
%     lincoeffs = 1./settings.orders;
    
    % estimate rho_on and gamma correction from second and third order
%     num = 50;
%     [gcor,rho] = kest_v2(sofi_c,num);
%     gammas = gcor./repmat((1:settings.orders(end))',1,size(gcor,2));
%     gammas(1,:)=1;
% 
% %     gammas = [1,1,1,1,1,1,1];
% 
%     figure, 
%     plot(mean(gammas,2));hold on;
%     plot(1./settings.orders,'r')
%     % gammas = 1./orders./0.9^(orders-1);
    
    fcpx = 0.65;
    iter = 100;
%     sofim=cellfun(@(x)squeeze(mean(reshape(x,size(x,1),size(x,2),5,[]),4)),sofi,'UniformOutput',0);
    lambda = [0.002,0.005,0.025,0.05];
    lambda = [0.002,0.01,0.025,0.05, 0.05, 0.05];
%     gammas = [0.7,0.4,0.3];
%     [sofi2]=sofiLinearize_mb2(sofi_c,settings.dec.fwhm,settings.sys.orders,iter,gammas,fcpx,1,lambda);
%     settings.dec.lambda = lambda;
    %Unload mex file - For testing purposes
%     clear mex;
    
    %% Save results
    
%     for order = settings.sys.orders
%         
%         temp = mean(sofi_c{order},3);
% %         temp = sum(sofi_c{order},3); % for PALM samples        
%         saveImageStack(temp,settings.io.outputpath,[settings.io.imageName,'sofi_raw',num2str(order)],[],16)
%         saveImageStack(sofi_c{order},settings.io.outputpath,[settings.io.imageName,'sofi_rawMovie',num2str(order)],[],16)
%         
%         temp = mean(sofi{order},3);
% %         temp = sum(sofi{order},3); % for PALM samples
%         saveImageStack(temp,settings.io.outputpath,[settings.io.imageName,'sofi_lucy',num2str(order)],[],16)
%         saveImageStack(sofi{order},settings.io.outputpath,[settings.io.imageName,'sofi_lucyMovie',num2str(order)],[],16)
%         
% %         temp = mean(sofif{order},3);
% %         temp = sum(sofi{order},3); % for PALM samples
% %         saveImageStack(temp,settings.io.outputpath,[settings.io.imageName,'sofi_lucyf',num2str(order)],[],16)
% %         saveImageStack(sofif{order},settings.io.outputpath,[settings.io.imageName,'sofi_lucyfMovie',num2str(order)],[],16)
%         
% %         temp = mean(sofi2{order},3);
% % %         temp = sum(sofi{order},3); % for PALM samples
% %         saveImageStack(temp,settings.io.outputpath,[settings.io.imageName,'sofi_breg',num2str(order)],[],16)
% %         saveImageStack(sofi2{order},settings.io.outputpath,[settings.io.imageName,'sofi_bregMovie',num2str(order)],[],16)
%         
% %         saveImageStack(sofi_lin{order},settings.io.outputpath,[settings.io.imageName,'sofi_lucy2Movie',num2str(order)],[],16)
% %         temp = mean(sofi_lin{order},3);
% %         temp = sum(sofi{order},3); % for PALM samples
% %         saveImageStack(temp,settings.io.outputpath,[settings.io.imageName,'sofi_lucy2',num2str(order)],[],16)
%         
%         
% %         temp = mean(sofi2{order},3);
% %         temp = sum(sofi2{order},3); % for PALM samples
% %         saveImageStack(sofi2{order},settings.outputpath,[settings.imageName,'sofi_bregMovie',num2str(order)],[],16)
% 
% %         temp = mean(sofi2{order},3);
% %         saveImageStack(temp,settings.outputpath,[settings.imageName, 'sofi_breg',num2str(order)],[],16)
%     end
    
    %% molecular parameters
    
    if settings.molpar.run ==1
        
%     if (~exist(outputPath,'dir'))
%     mkdir(outputPath);
%     end
        [ratio,density,brightness]=sofiParametersMy(sofi_c);
        sofimask = sofi_lin{4};
        alpha = mean(sofimask,3);
        alpha_max = max(alpha(:));
        alpha = imadjust(alpha./alpha_max,[settings.molpar.thresh 1],[0 1]);

        % alpha = alpha./max(alpha(:));
        alpha(alpha>0) = 1;
        density(isinf(density)) = 0;
    %     ratio = medfilt2(ratio,[3 3]);
        % ratio = imfilter(ratio,fspecial('gaussian',[3 3],1));

    %     brightness = medfilt2(brightness,[3 3]);
        % density = medfilt2(density,[2 2]);
        % density = imfilter(density,fspecial('gaussian',[9 9],1));

        I_ratio = ratio.*alpha;
        I_bright = brightness.*alpha;
        I_density = density.*alpha;

    %     figure, 
    %     subplot(311);imshow(I_ratio,[]);colormap('jet'); title('On time ratio');colorbar;
    %     subplot(312);imshow(I_bright,[]);colormap('jet'); title('Brightness');colorbar;
    %     subplot(313);imshow(I_density,[]);colormap('jet'); title('Number of molecules'); colorbar;

        saveImageStack(I_ratio,settings.io.outputpath,settings.io.imageName,'sofi_ratio',16);
        saveImageStack(I_bright,settings.io.outputpath,settings.io.imageName,'sofi_bright',16);
        saveImageStack(I_density,settings.io.outputpath,settings.io.imageName,'sofi_density',16);
%         save([settings.io.outputpath,filesep,settings.io.imageName,filesep,settings.io.imageName,'_molparams.mat'],'ratio','density','brightness','sofi_c','sofimask')
%         save([settings.io.outputpath,filesep,settings.io.imageName,filesep,settings.io.imageName,'_sofi_lin2.mat'],'sofi_lin')
%         save([settings.io.outputpath,filesep,settings.io.imageName,filesep,settings.io.imageName,'_sofi_lin1.mat'],'sofi')
    %     molStats(sofi_c,sofi,[settings.outputpath,filesep,settings.imageName],settings.wsize)
    end
end

%% Show results
% io = 2;

% figure,
% subplot(1,3,1);imshow(mean(stack,3),[]);title('Widefield');
% subplot(1,3,2);imshow(mean(sofi_c{io},3),[]);title('Cumulant image ');
% subplot(1,3,3);imshow(mean(sofi{io},3),[]);title('Lucy deconvolved dampar 0.0');
% fh = gcf;saveFigure(fh,outputpath,'wf_raw_lucy');

