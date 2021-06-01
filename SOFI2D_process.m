function [sofi,sofi_c,settings,sofi_lin,stats,results] = SOFI2D_process(settings, str) 
% 2D SOFI processing 
% Calculate SOFI cumulants, perform flattening, deconvolution and
% linearization according to bSOFI approach or adaptive linearization
% and additional computations like FRC resolution estimate and SNR estimation based on
% jacknife resampling and estimation of molecular parameters see bSOFI and PALM & SOFI paper
% Tomas Lukes, tomas.lukes@epfl.ch

    %
    %%% Load image stack
    %
    if strcmp(settings.io.fext,'.tif') == 1
        [stack,settings] = concatFiles(settings.io.imageName,settings,str);
        settings.sys.frames = size(stack,3);
           
    elseif strcmp(settings.io.fext,'.tiff') == 1
        [stack,frames] = loadStack(settings, str);
        settings.sys.frames = frames;   
    
    elseif strcmp(settings.io.fext,'.dat') == 1 || strcmp(settings.io.fext,'.raw') == 1 
        
        if strcmp(settings.io.fext,'.dat') == 1
            settings.io.byteorder = 'ieee-le'; 
        else
            settings.io.byteorder = 'ieee-be';
        end
        
        disp('Loading RAW data')
        
        filenamestr = char([settings.io.imageFile, str]);
        fname_temp = dir(filenamestr);
        %File lenght in frames based on its dimensions
        if isempty(settings.sys.end)== 0
            nframes = settings.sys.end;
        else
            nframes = fname_temp.bytes/(2*settings.io.W*settings.io.H);
        end
        
        disp(nframes);
        %nframes = settings.sys.sub;
        fid = fopen(filenamestr);
        fseek(fid, 0,'bof');
        stack = zeros(settings.io.W,settings.io.H,nframes, 'uint16');
        id = 1;
        for i=1:nframes
            data_temp = fread(fid,[settings.io.W,settings.io.H],'uint16',settings.io.byteorder);
            stack(:,:,id) = data_temp;
            id = id+1;
        end
        settings.sys.frames = nframes;
        settings.io.imageFile = char([settings.io.imagePath, filesep,settings.io.imageName]);
        settings.io.imageName = char(settings.io.imageName);
        fclose('all');
    end
    
  
    
    %%% Bleaching correction
    %
    [settings,results] = plotmitrace(stack,settings);

    if settings.io.blcor == 1
        [stack,results] = bleachcor(stack, settings,results);
    end
    
    %
    %%% Estimate Ton
    %
    if settings.ton.run ==1
        results = estimateton(stack,settings,results);
    end

    %% Drift correction 
    % drift correction based on tracking fiducial markers
    if settings.io.dcor == 1 
        if strcmp(settings.dcor.type,'TS')
            stack = driftcorTS(stack, settings,frames);
        elseif strcmp(settings.dcor.type,'LBEN_PALM')
            stack = driftcor(stack,settings);
        elseif strcmp(settings.dcor.type,'SOFI')
            disp('Drift correction from WF stack')
            fnumber = 500;
            wf1 = mean(stack(:,:,1:fnumber),3);
            
            ii_max = floor(settings.sys.frames/fnumber);
            drifts = zeros([settings.sys.frames, 2]);
            for ii = 1:ii_max
                wf = mean(stack(:, :, ((ii-1)*fnumber)+1:((ii-1)*fnumber)+fnumber),3);
                [drift_x, drift_y] = ccrShiftEstimation(wf1, wf, 2);
                if ii==ii_max
                    drifts(((ii-1)*fnumber)+1:end, :) = repmat([drift_x, drift_y],settings.sys.frames - (ii_max-1)*fnumber,1);
                else
                    drifts(((ii-1)*fnumber)+1:((ii-1)*fnumber)+fnumber, :) = repmat([drift_x, drift_y],fnumber,1);
                end
            end
            
            [sy,sx,~] = size(stack);
            [xx, yy]=meshgrid(1:sx,1:sy);

            for k=1:settings.sys.frames
                img = stack(:,:,k);

                imReg=interp2(xx,yy,single(img),xx-drifts(k,1),yy-drifts(k,2),'linear');
                imReg(isnan(imReg)) = 0;
                stack(:,:,k)=uint16(imReg);
            end
            
        else
            disp(['Wrong drift corr type :,',settings.dcor.type])
        end
    end

    %
    %%% Crop ROI
    %
    % crop roi after bleaching and drift correctons
    
    if isfield(settings.io,'roi')
		if ~isempty(settings.io.roi)
			roi = settings.io.roi;
			stack = stack(roi{2},roi{1},:);
		end
	end
    
    % remove frames if settings.sys.end is specified
        if ~isempty(settings.sys.end)
			stack = stack(:,:,1:settings.sys.end);
        end
    %end
	% remove frames if settings.sys.start is specified
	if ~isempty(settings.sys.start)
		stack = stack(:,:,settings.sys.start:end);
	end
		
    
    %% Calculate cumulants, perform flattening
    disp('sofiCumulantsAll')
    [sofi_c, settings,stats] = sofiCumulantsAll(stack, settings);

    % estmate drift from 2nd order sofi images
     if settings.io.dcor == 1 && strcmp(settings.dcor.type,'SOFI')
        im_fix = sofi_c{2}(:,:,1); 
        
        for ii = 1:size(sofi_c{2},3)
            
            % estimate the drift from 2nd order SOFI
            im_mov = sofi_c{2}(:,:,ii);
            [dx, dy] = ccrShiftEstimation(im_fix,im_mov,10);
            results.sofidrift(:,ii) = [dx, dy];%.*(settings.sys.pxy/2);
            
            if ii > 1 && settings.io.dcor == 1 
                % drift correction for all SOFI images
                for n = 2:numel(sofi_c)
                    [sy,sx,~] = size(sofi_c{n}(:,:,1));
                    [xx, yy]=meshgrid(1:sx,1:sy);
                    imReg=interp2(xx,yy,sofi_c{n}(:,:,ii),xx-n*dx/2,yy-n*dy/2,'linear');
                    imReg(isnan(imReg)) = 0;
                    sofi_c{n}(:,:,ii) = imReg;
                end
            end
        end
             
         figure, 
         plot(-results.sofidrift(1,:),'x-');hold on;
         plot(-results.sofidrift(2,:),'+-')
         ylabel('SOFI drift in pixels')
     end

    %
    %%% Linearization and postprocessing
    %
    disp('SOFIAdalin computation')
    sofi_lin = sofiAdalin(sofi_c, stack, settings); 
%     sofi_lin=0;

    if settings.dec.denoise == 1 
        [x,y]=xygrid(size(sofi_lin{2}(:,:,:)));
        temp=interpolate(sofi_c{2}(:,:,:),x,y);
    %     sofi_lin{2} = sum(temp.*sofi_lin{2}(:,:,:),3)./sum(temp,3);
        sofi_lin{2} = temp.*sofi_lin{2}(:,:,:)./repmat(sum(temp,3),1,1,size(sofi_lin{2},3));

        [x,y]=xygrid(size(sofi_lin{3}(:,:,:)));
        temp=interpolate(sofi_c{3}(:,:,:),x,y);
    %     sofi_lin{3} = sum(temp.*sofi_lin{3}(:,:,:),3)./sum(temp,3);
        sofi_lin{3} = temp.*sofi_lin{3}(:,:,:)./repmat(sum(temp,3),1,1,size(sofi_lin{3},3));

        [x,y]=xygrid(size(sofi_lin{4}(:,:,:)));
        temp=interpolate(sofi_c{4}(:,:,:),x,y);
    %     sofi_lin{4} = sum(temp.*sofi_lin{4}(:,:,:),3)./sum(temp,3);
        sofi_lin{4} = temp.*sofi_lin{4}(:,:,:)./repmat(sum(temp,3),1,1,size(sofi_lin{4},3));

    end
    
    disp('sofiLinearize')
    sofi=sofiLinearizeCooseDeconvolution(sofi_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,settings);
%     sofi = 0;
    
    
    %% Save SOFI results

%     sofiSaveResults(settings,sofi_c,sofi)
    sofiSaveResults(settings,sofi_c, sofi_lin, sofi)
    
    %
   %% Estimate molecular parameters
    
    if settings.molpar.run ==1 && sum(settings.sys.orders) >= 9 
        mparams = molparams(sofi_c,sofi,settings);
        results.mparams = mparams;
    end
    
    if isfield(settings,'rep') && settings.rep.run == 1
        imIn = mean(stack,3);
        printReport(settings,results,sofi_c,sofi_lin,imIn,mparams);
    end
    
    %
    %% FRC estimation
    %
    if settings.frc.run == 1
        estimResol(sofi_lin, settings);
    end

%% Show results
% io = 2;

% figure,
% subplot(1,3,1);imshow(mean(stack,3),[]);title('Widefield');
% subplot(1,3,2);imshow(mean(sofi_c{io},3),[]);title('Cumulant image ');
% subplot(1,3,3);imshow(mean(sofi{io},3),[]);title('Lucy deconvolved dampar 0.0');
% fh = gcf;saveFigure(fh,outputpath,'wf_raw_lucy');

