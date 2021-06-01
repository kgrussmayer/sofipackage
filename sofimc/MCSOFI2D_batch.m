%%% batch code for MCSOFI processing
%%% adrien.descloux@epfl.ch 11/2018
%%% tomas.lukes@epfl.ch
%%% kristin.grussmayer@epfl.ch

close all; clear; 

%
%%% Initialization

addpath(genpath('..\utils'));
addpath(genpath('funcs'));
config_MCSOFI2D;

% select folder to batch process

% pname = uigetdir('K:\Kristin');
% name of the folder that contains the data files and the calibration
% folder
pname = 'K:\Kristin\20180201-selectionLagTime1'; 
settings.io.pn = pname;
% name of the folder under pname that contains the calibration files
calibrationFolderName = 'calibration-BS640';

list = dir(pname);
list(1:2) = []; % remove . and .. path from the list

% make new results folder
rFold = ['results_',getID];
mkdir(pname,rFold);

% look for calibration file
isCal = 0;
for k = length(list):-1:1 % loop backwards because we are changing the list on the fly
    if strfind(list(k).name,calibrationFolderName)
        settings.io.pnc = [pname,filesep,list(k).name];
        isCal = 1;
        list(k) = []; % remove calibration folder from the list
    elseif strfind(list(k).name,'results')
        list(k) = []; % remove results folder from the list
    elseif isempty(strfind(list(k).name,'.tiff'))&& isempty(strfind(list(k).name,'.tif'))
        list(k) = []; % remove extra folders from tomas procesing
    elseif strfind(list(k).name,'.mat')
        list(k) = [];
    end
end
if isCal == 0
    disp('Calibration file not found !');
    disp('Please make sure the folder is properly organized!'); 
else

% make list of paired files
ex = length(list(1).name)-strfind(list(1).name,'.');
idMax = 0; listC = [];
for k = 1:length(list)
    temp = str2double(list(k).name(end-ex-3:end-ex-1));
    if temp > idMax
        idMax = temp;
    end
end

rList = []; tList = [];
for k = 0:idMax 
    id = num2str(k); while length(id) < 3; id = ['0',id];end
    id = ['_',id];
    for h = length(list):-1:1
       if strfind(list(h).name,id) % if we find proper ID number
           if strfind(list(h).name,'reflected')
               rList(end+1).name = list(h).name;
               rList(end).id = id(2:end);
           elseif strfind(list(h).name,'transmitted')
               tList(end+1).name = list(h).name;
               tList(end).id = id(2:end);
           end
       end
    end
end

if length(tList) ~= length(rList)
    warning('Error number of reflected and transmitted files are different');
end

% summary of folder to be processed
disp(['Selected folder : ',settings.io.pn])
disp(['Calibration folder : ',settings.io.pnc])
disp(['Number of files detected : ',num2str(length(rList))])
end

%% perform beads based calibration if selected in settings and if calibration files were found

if settings.cal.beads == 1  
    if isCal == 0
        disp('Calibration file not found !');
        disp('Please make sure the folder is properly organized!'); 
    else
        %%% Load the calibration files
        imc1 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc1, settings.io.fileExtension]);
        imc2 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc2, settings.io.fileExtension]);
        if strcmp(settings.io.dataType, 'exp')
            imc2 = flip(imc2,2);
        end

        % imc1(imc1==0)=min(imc1(imc1>0));
        % imc2(imc2==0)=min(imc2(imc2>0));

        % what is indcom used for?
        [~,indcom]=max(squeeze(mean(max(imc1))).*squeeze(mean(max(imc2))));
        im_fix=imc1(:,:,indcom);
        im_mov=imc2(:,:,indcom);

        show2subs(im_fix,im_mov,1)
        close(1)

        %%% Run calibration (register images from 2 color channels)
        disp('Compute calibration')
        settings.cal.roix = []; settings.cal.roiy = [];
        cal = runcalib_mcsofi(imc1,imc2,settings.sys,settings.cal,settings.io); 
    end
end

%% MAIN PROCESSING LOOP
for k = 1:length(rList)
    disp(['Processing file # ',num2str(k),', ID : ',rList(k).id])
    %% Load stack of images to be transformed (1st stack) 
    settings.io.fn1 = rList(k).name;
    settings.io.fn2 = tList(k).name;
    disp('Loading reflected data')
    imstack1 = load_tifFile([settings.io.pn, filesep, settings.io.fn1],settings.sys.sub);
    disp('Loading transmitted data')
    imstack2 = load_tifFile([settings.io.pn, filesep, settings.io.fn2],settings.sys.sub);
    if strcmp(settings.io.dataType, 'exp')
        imstack2=flip(imstack2,2);
    end

    disp('Compute stacks STD')
    % compute the std of the stacks, should change to start frames : end
    % frames
    st1 = std(single(imstack1(:,:,1000:end-1000)),[],3);
    st2 = std(single(imstack2(:,:,1000:end-1000)),[],3);

    %% use cross-correlations between std of acquired image sequences to
    % calculate the transformation
    if settings.cal.beads == 0
        disp('Compute correct coregistration')
        % 2D high pass filter
        kernel1 = ones(3)/9;
        t1 = st1(50:end-50,50:end-50); t2 = st2(50:end-50,50:end-50); 
        t1 = medfilt2(t1,[2 2]); t2 = medfilt2(t2,[2 2]); % remove salt&pepper noise
        t1 = t1-imgaussfilt(single(t1),5); t2 = t2-imgaussfilt(single(t2),5); % remove background before cross-correlation

        if strcmp(settings.io.dataType, 'exp')

            temp = xcorr2(t1,t2);
            temp((size(temp,1)+1)/2,:) = 0;
            temp(:,(size(temp,2)+1)/2) = 0;
            [~,ind] = max(temp(:));
            [y,x] = ind2sub(size(temp),ind);
            shiftx = -(x-(size(temp,2)+1)/2);
            shifty = -(y-(size(temp,1)+1)/2);
        else
            shiftx = 0;
            shifty = 0;
        end
    end
    %% Transform stack of images with the chosen calibration method
    disp('Transform of imstack2')
    Nframes = size(imstack1,3);
    imstack2t = zeros(size(imstack1),'single');
    % apply transform to data
    if settings.cal.beads == 0 % using std of raw data
        for ii = 1:Nframes

            im_mov = imstack2(:,:,ii);    
            im1t = imtranslate(im_mov,[-shiftx -shifty]);
            imstack2t(:,:,ii) = single(im1t);
        end
        st2t = imtranslate(st2,[-shiftx -shifty]);
        settings.cal.shiftx = -shiftx;
        settings.cal.shifty = -shifty; % 
    else % using beads calibration
        for ii = 1:Nframes
            im_mov = imstack2(:,:,ii);   
            im1t = imwarp(im_mov,cal.tf{1},'OutputView',imref2d(size(im_mov)));
            imstack2t(:,:,ii) = single(im1t);
        end
        st2t = imwarp(st2,cal.tf{1},'OutputView',imref2d(size(st2)));
        settings.cal.shiftx = -cal.tf{1}.A(1);
        settings.cal.shifty = -cal.tf{1}.B(1);
    end

    % compute coregistration mask and proper crop
    mask = not(im1t(:,:,1) == 0);
    [r,c] = find(mask>0);
    rect = [min(c)+1 min(r)+1 max(c)-min(c)-2 max(r)-min(r)-2];
    % SOFI processing requires even dimensions or you gets weird cropping
    if mod(rect(3),2) == 0; rect(3) = rect(3)-1; end
    if mod(rect(4),2) == 0; rect(4) = rect(4)-1; end

    settings.cal.roix = [rect(1) rect(3)];
    settings.cal.roiy = [rect(2) rect(4)];

    % crop image stack 
    imstack1crop = imstack1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);
    % crop image stack 
    imstack2crop = imstack2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);

    % calculate the mean transformed image
    mean1c = mean(double(imstack1crop), 3);
    mean2tc = mean(double(imstack2crop), 3);
    % crop std 
    st1c = st1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3));
    st2tc = st2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3)); 
    fake = [];
    fake(:,:,2) = 255.*st1c./max(st1c(:)); 
    fake(:,:,3) = 255.*st2tc./max(st2tc(:)); 
    fake(:,:,1) = 255.*st2tc./max(st2tc(:));
    figure(2);imagesc(uint8(fake));title('Std overlay');
    if settings.io.figshow == 0
        close
    end
        
    %% Calculate MC SOFI cumulants of cropped raw data (1st stack) and transformed and cropped 2nd stack
    disp('Compute SOFI cumulants')
    orders = settings.sys.orders;
    subseqlength = settings.sys.wsize;
    start = settings.sys.start;

    Nss=floor((Nframes-start+1)/subseqlength); % number of subsequences

    for ns=1:Nss
    %     disp(['processing Nss:',num2str(ns)]);
        fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
        sofi=mcSofiCumulants(permute(cat(4,imstack1crop(:,:,fr),imstack2crop(:,:,fr)),[1 2 4 3]),[],[],[],orders);
        if ns==1
            c=sofi;
            c=cellfun(@(x)repmat(x,[1 1 1 Nss]),c,'UniformOutput',0);
        else
            for io=orders
                c{io}(:,:,:,ns)=sofi{io};
            end
        end
    end

    %% Calculate 2D cross-cumulants of non-transformed raw data (2nd stack)
    disp('Compute SOFI cumulant of imstack2')
    for ns=1:Nss
    %     disp(['processing Nss:',num2str(ns)]);
        fr=start+(ns-1)*subseqlength-1+(1:subseqlength);
        sofi=mcSofiCumulantsSingle(imstack2(:,:,fr),[],[],[],orders);
        if ns==1
            c2=sofi;
            c2=cellfun(@(x)repmat(x,[1 1 Nss]),c2,'UniformOutput',0);
        else
            for io=orders
                c2{io}(:,:,ns)=sofi{io};
            end
        end
    end

    %% Flattening before unmixing if chosen
    if strcmp(settings.sys.flattening, 'before')
        disp('Cumulant flattening')
        c = sofiAllFlatten3D(c, orders);

        % flattening of 2D cross-cumulants of non-transformed raw data
        c2 = sofiAllFlatten(c2,orders);
    end

    %% Transform second channel cross-cumulants and combine cumulants
    disp('Transform and crop second SOFI channel')
    c2t=c2;
    off = [1 4 6];
    for io=orders
    %     disp(io)
    %     T=tform.tdata.T;
    %     T(3,1:2)=io*T(3,1:2);
    %     T(3,2)=T(3,2)-io; %might have to be corrected differently - case by case
    %     c2t{io}=imtransform(c2{io},maketform('affine',T),'XData',[1 size(c{io},2)],'YData',[1 size(c{io},1)],'Size',[size(c{io},1) size(c{io},2)]);
    %      c2t{io} = imwarp(c2{io}, cal.tf{io},'OutputView',imref2d(size(c2t{io})));
        if settings.cal.beads == 0 
            % transform
            c2t{io} = imtranslate(c2{io},io.*[settings.cal.shiftx settings.cal.shifty]); % -cal.tf{1}.B(1)
            % crop
            c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-off(io):io*rect(2) + io*rect(4)-ceil(3*io/2)-off(io),...
                 io*rect(1)+floor(3*io/2)-off(io):io*rect(1)+io*rect(3)-ceil(3*io/2)-off(io),:);
        else% using beads calibration   
            % transform % check 
            c2t{io} = imwarp(c2{io},cal.tf{io},'OutputView',imref2d(size(c2{io})));
            % crop
            c2t{io} = c2t{io}(io*rect(2)+floor(3*io/2)-off(io):io*rect(2) + io*rect(4)-ceil(3*io/2)-off(io),...
                 io*rect(1)+floor(3*io/2)-off(io):io*rect(1)+io*rect(3)-ceil(3*io/2)-off(io),:);
        end
    end

    disp('Combine cumulants')
    tmp=cellfun(@(x)reshape(x(:,:,:),size(x,1),size(x,2),1,[]),c2t,'UniformOutput',0);
    for io=orders
        c{io}(:,:,end,:)=tmp{io};
    end

    %% Merging second order cross-cumulants into RGB image
    % choose second order
    mc2=c{2};
    im_rgb_cc2 = mergeToRgb(mean(squeeze(mc2(:,:,1,:)),3), mean(squeeze(mc2(:,:,2,:)),3), mean(squeeze(mc2(:,:,3,:)),3));

    figure(3);
    image(im_rgb_cc2);
    title('Second order cross-cumulants');
    set(gca,'xtick',[],'ytick',[]);
    if settings.io.figshow == 0
        close
    end
    %% Spectral unmixing - using 2nd order cumulant - 3 colors
    disp('Spectral unmixing')
    % load the transmission coefficient T (proportion in transmission channel)
    T1 = settings.mc.T1;
    T2 = settings.mc.T2;
    T3 = settings.mc.T3;
    % load the ratio R, this is not the reflection coefficient!!
%     R1 = settings.mc.R1;
%     R2 = settings.mc.R2;
%     R3 = settings.mc.R3;

    if settings.mc.getUnmixCoefAlexa488XAlexa647 == 1
        % option to read coefficients from file Alexa488XAlexa647!
        % get T1,T2,T3 from calibration and filename information
        if strcmp(settings.io.dataType,'exp')
            [T1, T2, T3] = getUnmixCoefAlexa488XAlexa647Exp(settings.io.pnc,settings.io.fn1);
        else
            [T1, T2, T3] = getUnmixCoefAlexa488XAlexa647Theo(settings.io.pnc,settings.io.fn1);
        end
        settings.mc.T1 = T1;
        settings.mc.T2 = T2;
        settings.mc.T3 = T3;
        disp('Warning: loads unmixing coefficients from file. Assuming combination Alexa488 X Alexa647!')
    end

%     R=[R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1]; % old formulation
    R = [T1^2 T2^2 T3^2; 
        T1*(1-T1) T2*(1-T2) T3*(1-T3); 
        (1-T1)^2 (1-T2)^2 (1-T3)^2];
    Rinv=R^(-1);

    mc2=c{2};
    umx1c2=squeeze(Rinv(1,1)*mc2(:,:,3,:)+Rinv(1,2)*mc2(:,:,2,:)+Rinv(1,3)*mc2(:,:,1,:));
    umx2c2=squeeze(Rinv(2,1)*mc2(:,:,3,:)+Rinv(2,2)*mc2(:,:,2,:)+Rinv(2,3)*mc2(:,:,1,:));
    umx3c2=squeeze(Rinv(3,1)*mc2(:,:,3,:)+Rinv(3,2)*mc2(:,:,2,:)+Rinv(3,3)*mc2(:,:,1,:));
    %% Flattening after unmixing if chosen
    if strcmp(settings.sys.flattening, 'after')
        % prepare to also save the unflattened, unmixed second order data
        im_rgb_sofi2NoFlat = mergeToRgb(mean(umx3c2,3), mean(umx2c2,3), mean(umx1c2,3));
        % flattening after unmixing
        umx_c2 = {umx1c2, umx2c2, umx3c2};
        umx_c2_flat = sofiAllFlattenMC(umx_c2,2);
        umx1c2 = umx_c2_flat{1};
        umx2c2 = umx_c2_flat{2};
        umx3c2 = umx_c2_flat{3};
    end

    %% Merging unmixed cumulants into RGB image
    im_rgb_sofi2 = mergeToRgb(mean(umx3c2,3), mean(umx2c2,3), mean(umx1c2,3));

    figure(4);
    image(im_rgb_sofi2);
    title('Unmixed second order cumulants');
    set(gca,'xtick',[],'ytick',[]);
    if settings.io.figshow == 0
        close
    end

    %% unmix 3 colors 3rd order

    % R=[R1^3 R2^3 R3^3; R1^2 R2^2 R3^2; R1 R2 R3; 1 1 1]; old formulation
    R = [   T1^3        T2^3      T3^3      ;
         T1^2*(1-T1) T2^2*(1-T2) T3^2*(1-T3);
         T1*(1-T1)^2 T2*(1-T2)^2 T3*(1-T3)^2;
          (1-T1)^3     (1-T2)^3   (1-T3)^3];

    Rinv=pinv(R);

    mc3 = c{3};
    nth=zeros(size(mc3,1),size(mc3,2),2);
    umx1c3=squeeze(Rinv(1,1)*mc3(:,:,4,:)+Rinv(1,2)*mc3(:,:,3,:)+Rinv(1,3)*mc3(:,:,2,:)+Rinv(1,4)*mc3(:,:,1,:));
    umx2c3=squeeze(Rinv(2,1)*mc3(:,:,4,:)+Rinv(2,2)*mc3(:,:,3,:)+Rinv(2,3)*mc3(:,:,2,:)+Rinv(2,4)*mc3(:,:,1,:));
    umx3c3=squeeze(Rinv(3,1)*mc3(:,:,4,:)+Rinv(3,2)*mc3(:,:,3,:)+Rinv(3,3)*mc3(:,:,2,:)+Rinv(3,4)*mc3(:,:,1,:));

    %% flattening after unmixing
    if strcmp(settings.sys.flattening, 'after')
        umx_c3 = {umx1c3, umx2c3, umx3c3};
        umx_c3_flat = sofiAllFlattenMC(umx_c3,3);
        umx1c3 = umx_c3_flat{1};
        umx2c3 = umx_c3_flat{2};
        umx3c3 = umx_c3_flat{3};
    end

    im_rgb_sofi3 = mergeToRgb(umx3c3, umx2c3, umx1c3);
    im_rgb_sofi3(im_rgb_sofi3>1)=1;

    %% Deconvolution and Linearization of second order unmixed and flattened SOFI images
    disp('Deconvolution')
    [sy,sx,~] = size(umx1c2);
    order=2;

    % matrix rows, cols, frames, color channels
    umx_c2_all = cat(4,umx1c2, umx2c2, umx3c2);
    im_out = [];

    % h2 = fspecial('Gaussian',[round(sy/2),round(sx/2)],2.3*sqrt(order)/2/sqrt(2*log(2))*sy);
    h2 = fspecial('Gaussian',[sy,sx],settings.dec.fwhm*sqrt(order)/sqrt(8*log(2)));
    % h2 = padarray(h2,[sy-round(sy/2),sx-round(sx/2)]);
    h = fspecial('Gaussian',29,settings.dec.fwhm*sqrt(order)/sqrt(8*log(2)));

    %DEFINE DECONV PARAMS
    % parameters for cuda version of bregman iterative method based
    % deconvolution with Gaussian noise model
    settings.iter = 50;
    settings.apodize =0;
    settings.lambda = 0.03;
    settings.omega = 0.95;
    settings.NumImages = 1;

    % parameters for matlab version of augmented lagrangian based deconvolution
    % with Gaussian noise model
    settings.gamma = 2000;
    settings.beta = settings.gamma;
    settings.alpha = 1;
    settings.reltol = 1e-4;
    settings.maxIter = 20;
    settings.Lp = 1;

    for n=1:size(umx1c2,3)
        for index_colors = 1:size(umx_c2_all,4)
            tmp = umx_c2_all(:,:,n,index_colors);
            tmp(tmp<0) = 0;
            scale=max(tmp(:));
            tmp = tmp./scale;
            % deconvolution depending on the chosen method
            if strcmp(settings.deconv_method,'breg_cuda')
                %%% CUDA BEGIN S2D %%%
                Aaray = [settings.iter,settings.apodize,settings.lambda,...
                    settings.omega,settings.NumImages];
                im_dec = SOFI_X7(tmp,h2,Aaray);
                %%% CUDA END S2D %%%
            elseif strcmp(settings.deconv_method,'augLag_mb')

                [imOut,relChangeAll] = deconvAugLag2D(tmp,fftn(h,[sy,sx]),settings);
                im_dec = imOut;
            elseif strcmp(settings.deconv_method,'lucy')
                im_dec = deconvlucy(tmp,h,settings.dec.iter);
                
            elseif strcmp(settings.deconv_method,'bsofi')
                im_dec = deconvlucy(tmp,h,50);
            else
                disp('Deconvolution method not implemented')
            end
            % optional linearization
            if settings.dec.lin == 1;
                im_out(:,:,n,index_colors) = scale*(im_dec.^settings.sys.lincoeff(2));
            else
                im_out(:,:,n,index_colors) = scale*(im_dec);
            end
            
            if strcmp(settings.deconv_method,'bsofi');
                % reconvolve the image
                psf=gaussian(settings.dec.fwhm/sqrt(8*log(2))); % 
                psf=psf(:)*gaussian(settings.dec.fwhm/sqrt(8*log(2)));
                im_out(:,:,n,index_colors)=convn(im_out(:,:,n,index_colors),psf,'same');   
            end
        end
    end

    %% Merging linearized deconvolved unmixed cumulants into final RGB image
    im1 = mean(im_out(:,:,:,1),3);
    im2 = mean(im_out(:,:,:,2),3);
    im3 = mean(im_out(:,:,:,3),3);
    im_rgb_sofi2_lin = mergeToRgb(im3, im2, im1);

    r_lin = im_rgb_sofi2_lin(:,:,1);
    g_lin = im_rgb_sofi2_lin(:,:,2);
    b_lin = im_rgb_sofi2_lin(:,:,3);

    figure(5)
    subplot(131);
    imshowpair(mean1c,mean2tc,'falsecolor');axis equal;axis tight;
    title('Mean overlay')
    subplot(132);
    imshow(im_rgb_sofi2,[])
    title('SOFI 2 unmixed')
    subplot(133);
    imshow(im_rgb_sofi2_lin,[]);
    title('SOFI 2 lin dec unmixed')
    set(gcf,'position',[465   555   1200   350]);
    if settings.io.figshow == 0
        set(gcf,'visible', 'off');
    end

    figure(6)
    subplot(131);
    rf = zeros(size(r_lin,1),size(r_lin,2),3) ; rf(:,:,1) = r_lin;  rf(:,:,3) = r_lin; 
    imshow(rf./max(rf(:)));
    title('Red channel')
    subplot(132);
    gf = zeros(size(g_lin,1),size(g_lin,2),3) ; gf(:,:,2) = g_lin;  gf(:,:,1) = g_lin; 
    imshow(gf./max(gf(:)));
    title('Green channel')
    subplot(133);
    bf = zeros(size(b_lin,1),size(b_lin,2),3) ; bf(:,:,3) = b_lin;  bf(:,:,2) = b_lin; 
    imshow(bf./max(bf(:)));
    title('Blue channel')
    set(gcf,'position',[465   555   1200   350]);
    if settings.io.figshow == 0
        set(gcf,'visible', 'off');
    end

    %% saving 
    if settings.io.figsave ==1 
        disp('Saving results')
        settings.io.outputpath = [pname,filesep,rFold];

        % save std and mean of stacks 
        writeTIFF(st1c./max(st1c(:)),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_std_R'])
        writeTIFF(st2tc./max(st2tc(:)),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_std_T'])
        writeTIFF(mean1c./max(mean1c(:)),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_mean_R'])
        writeTIFF(mean2tc./max(mean2tc(:)),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_mean_T'])
        writeRGBTIFF(uint8(fake),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_coreg'])

        % should include saving of cross-cumulants before unmixing
        writeTIFF(im_rgb_cc2,[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_cc2']);
        writeRGBTIFF(uint8(255.*im_rgb_cc2),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_cc2_overlay'])

        if strcmp(settings.sys.flattening, 'after')
            % save unmixed 2nd before flattening
            writeTIFF(im_rgb_sofi2NoFlat,[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixedNoFlat2']);
            writeRGBTIFF(uint8(255.*im_rgb_sofi2NoFlat),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixedNoFlat2_overlay'])
        end

        % save unmixed 2nd
        writeTIFF(im_rgb_sofi2,[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixed2']); 
        writeRGBTIFF(uint8(255.*im_rgb_sofi2),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixed2_overlay'])

        % save unmixed 3rd
        writeTIFF(im_rgb_sofi3,[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixed3']); 
        writeRGBTIFF(uint8(255.*im_rgb_sofi3),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixed3_overlay'])

        % save lin deconv
        writeTIFF(im_rgb_sofi2_lin,[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixed2_dec_lin']); 
        writeRGBTIFF(uint8(255.*im_rgb_sofi2_lin),[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_unmixed2_dec_lin_overlay'])

        saveFigure(5,[settings.io.pn,filesep,rFold],[rList(k).id,'_results'],'tif')
        saveFigure(6,[settings.io.pn,filesep,rFold],[rList(k).id,'_results_RGB'],'tif')
        saveSettingsTxt(settings,[settings.io.pn,filesep,rFold,filesep,rList(k).id,'_settings.txt'])
    end

end




