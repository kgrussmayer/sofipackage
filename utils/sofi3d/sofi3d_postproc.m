function [sofizt,sofid,dec] = SOFI3D_postproc(sofizt,dec,cal,io)
% Multiplane SOFI - processing raw 3D cumulants
% author: Tomas Lukes, tomas.lukes@epfl.ch
% MATLAB version 2012a

% note: This script run the processing of raw 3D cumulants 
% 1) Flattening and registration of 3D cumulants
% 2) Deconvolution and linearization
% 3) Saving the results
 

    
%% Perform flattening 

    sofizt = sofiAllFlatten3D(sofizt,dec.orders);
    
%% Register 3D cumulants

%     tfmat = out.tf; % !!!
    tfmat = cal.tf;  % Transformation matrices
    disp('Register 3D cumulants')
    for order=dec.orders
        for m=0:dec.nplanes-3  
            T=1;
            for mi=0:m
                T=T*tfmat{mi+1}.tdata.T;
            end
    %         T=T*[order 0 0; 0 order 0; 0 0 1];
            T(3,1:2)=T(3,1:2)*order;%/3;
%             tftmp=maketform('affine',T);

            tmp=sofizt{order}(:,:,(m+1)*order+(1:order),:);
%             tmp=imtransform(tmp, tftmp, 'XData', [1 size(tmp,2)], 'YData', [1 size(tmp,1)], ...
% %                'Size', [size(tmp,1),size(tmp,2)]);
            tmp=imwarp(tmp, imref2d(size(squeeze(tmp))),affine2d(T),...
                'OutputView',imref2d(size(squeeze(tmp))));
            sofizt{order}(:,:,(m+1)*order+(1:order),:)=tmp;
        end
    end
    
    % coregistration of the last plane
    for order=dec.orders
        m=dec.nplanes-2;
        T=1;
        for mi=0:m
            T=T*tfmat{mi+1}.tdata.T;
        end
    %         T=T*[order 0 0; 0 order 0; 0 0 1];
        T(3,1:2)=T(3,1:2)*order;%/3;
%         tftmp=maketform('affine',T);

        tmp=sofizt{order}(:,:,(m+1)*order+1,:);
%         tmp=imtransform(tmp, tftmp, 'XData', [1 size(tmp,2)], 'YData', [1 size(tmp,1)], ...
%            'Size', [size(tmp,1),size(tmp,2)]);
        tmp=imwarp(tmp, imref2d(size(squeeze(tmp))),affine2d(T),...
                'OutputView',imref2d(size(squeeze(tmp))));
        sofizt{order}(:,:,(m+1)*order+1,:)=tmp;
    end
    
%% Crop data if settings.sys.cropData == 1
if dec.cropData == 1 % crop data to keep only usable 3D data based on beads coregistration
    disp('Crop data based on coregistration')
    for order=dec.orders
        dx = order*cal.dx; Lx = order*cal.Lx;
        dy = order*cal.dy; Ly = order*cal.Ly;
        if size(sofizt{order},1) > dy+Ly
            sofizt{order} = sofizt{order}(dy:dy+Ly,dx:dx+Lx,:,:);
        end
    end
end

%% post processing 

if dec.driftCorr == 1
    disp('Drift correction based on sofi2 data')
    % estimate drift based on sofi2 data
    shiftx = []; shifty = [];
    sofi2 = squeeze(sofizt{2}(:,:,4,:));

    for n = 1:size(sofizt{2},3)
        sofi2 = squeeze(sofizt{2}(:,:,n,:));
        im_fix = sofi2(:,:,1);
        for k = 2:size(sofi2,3)
            im_mov = sofi2(:,:,k);
            [shiftx(n,k),shifty(n,k)] = ccrShiftEstimation(im_fix,im_mov,4);
        end
    end
    shiftx = median(shiftx,1);
    shifty = median(shifty,1);
    disp(['x drift estimation: ',num2str(shiftx,3), ' (sofi2 pixels)'])
    disp(['y drift estimation: ',num2str(shifty,3), ' (sofi2 pixels)'])
    % apply the drift correction
    for order = dec.orders
        for k = 2:size(sofi2,3)
            for n = 1:size(sofizt{order},3)
                sofizt{order}(:,:,n,k) = imtranslate(sofizt{order}(:,:,n,k),[order*shiftx(k)/2,order*shifty(k)/2]);
            end
        end
    end
end

if dec.avgn >= 1
    disp('Sofi sequence averaging')
% sofizt=cellfun(@(x)squeeze(mean(reshape(x,size(x,1),size(x,2),size(x,3),dec.avgn,[]),4)),sofizt,'UniformOutput',0);
    for order = dec.orders
    	sofizt{order} = mean(sofizt{order},4);
    end
end

if dec.medfilt == 1
    disp('2D median filtering')
    % 2D median filtering
    for o = 1:numel(sofizt)
        for n=1:size(sofizt{o},3)
            for m=1:size(sofizt{o},4)
                sofizt{o}(:,:,n,m)=medfilt2(sofizt{o}(:,:,n,m),[3 3]);
            end
        end
    end
end

%% Resolution estimate
if dec.estimRes == 1 && exist('../ImDecorr','dir') >= 1
    disp('Estimate the resolution')
    addpath(genpath('../ImDecorr')); % put the whole toolbox in Matlab path
    dec.resx = [];
    dec.resz = []; % there are no tools for axial resolution estimate yet
    for order = dec.orders
        % average the stack in a single frame
    	temp = mean(mean(sofizt{order},3),4);
        temp = apodImRect(temp(1:min(size(temp)),1:min(size(temp))),20);
        dec.resx(order) = 2/getDcorr(gpuArray(temp));
        disp(['Estimated resolution, order ',num2str(order),' ,',num2str(dec.resx(order),3),' (pixel)'])
    end
end
if exist('../ImDecorr','dir') == 0
    disp('Resolution estimate toolbox not on the expected path "../ImDecorr"')
end

%% Deconvolution
disp('Deconvolution')
[sofid]=sofiLinearize3Dpad(sofizt,dec);

% dec.lincoeffs = [1 0.5 0.33 0.25];
% dec.fcpx = 0.65;
% dec.iter = 50;
% [sofid]=sofiLinearize_mb3D(sofizt,dec);

% ADD OTHERS DECONVOLUTION AugLag/ Bregman

%% Save results
% ii = 2;

if (~exist(io.outputpath,'dir'))
    mkdir(io.outputpath);
    mkdir([io.outputpath]);
end

outputpath = [io.outputpath,filesep,'results'];

for ii = dec.orders
    
    sofi_l = mean(sofid{ii},4);
    saveImageStack(sofi_l,outputpath,[io.fileName,'_ord_',num2str(ii),'_lucy'],[],16);

    for kk = 1:io.seq
        sofi_l = sofid{ii}(:,:,:,kk);
        saveImageStack(sofi_l,outputpath,[io.fileName,'_ord_',num2str(ii),'_lucy_f',num2str(kk)],[],16);

        sofi_l = sofizt{1}(:,:,:,kk);
        saveImageStack(sofi_l,outputpath,[io.fileName,'_ord_',num2str(ii),'raw_f',num2str(kk)],[],16);
    end
    
    sofi_c = mean(sofizt{ii},4);
    saveImageStack(sofi_c,outputpath,[io.fileName,'_ord_',num2str(ii),'_sofic16'],[],16);
end     
