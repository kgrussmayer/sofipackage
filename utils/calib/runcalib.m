function cal = runcalib(sys,cal,pnc)
    cfext = sys.fext;
    tf = [];
    
    if strcmp(cfext,'.bin')
        file1='data1.bin';
        file2='data2.bin';

        load([pnc filesep 'settings.mat']);

        fid1=fopen([pnc filesep file1]);
        fid2=fopen([pnc filesep file2]);
        cfinfo=dir([pnc filesep file1]);

        nx=cam1.ROIPosition(3);
        ny=cam1.ROIPosition(4);

        N=cfinfo.bytes/2/nx/ny;

        k=0;
        data1c=uint16(zeros(ny,nx,N));
        data2c=data1c;
        while k<N-1
            k=k+1;
            data1c(:,:,k)=fread(fid1,[ny,nx],'*uint16');
            data2c(:,:,k)=fread(fid2,[ny,nx],'*uint16');
        end
        clear k;
        fclose(fid1);
        fclose(fid2);
        
    elseif strcmp(cfext,'.dat')
        calibration_file = [pnc,filesep,'beads'];
        % try calibration dat file
        if exist([calibration_file, '.dat'], 'file') 
            fidc = fopen([calibration_file, '.dat']);
            cfinfo=dir([calibration_file, '.dat']);
            N = cfinfo.bytes /(2*2*600*2048);
            data1c=uint16(zeros(600,2048,N));
            data2c=data1c;
            id = 1;
            for k = 1:(2*N)
                if mod(k,2) == 1
                    temp = fread(fidc,[2048,600],'uint16','ieee-be');
                    data1c(:,:,id) = temp';
                else
                    temp = fliplr(fread(fidc,[2048,600],'uint16','ieee-be'));
                    data2c(:,:,id) = rot90(temp',2);
                    id = id +1; 
                end
            end
            fclose(fidc);
            
        % try calibration tiff file
        elseif exist([calibration_file, '.tif'], 'file')
            
            dataTemp = load_tiff([calibration_file, '.tif']);
            data1c = dataTemp(:,:,2:2:end);
            data2c = dataTemp(:,end:-1:1,1:2:end-1);
        else
            errID = 'MATLAB:inputError';
            msgtext = 'Calibration file %s does not exist';
            ME = MException(errID,msgtext, [calibration_file, '.tif']);
            throw(ME)
        end
            
        
    elseif strcmp(cfext,'.tif') || strcmp(cfext,'.tiff')

        dataTemp = load_tiff([pnc,filesep,'beads.tif']);
        data1c = dataTemp(:,:,2:2:end);
        data2c = dataTemp(:,end:-1:1,1:2:end-1);
    end

    %%% Reorder planes for data using new prism
    if strcmp(cfext,'.bin')
        data1c = reorderData(data1c,1,0,sys);
        data2c = reorderData(data2c,2,1,sys);

        [sy,sx,st] = size(data1c);
        data1c = reshape(data1c,sy,sx/4,4,st);
        data2c = reshape(data2c,sy,sx/4,4,st);
        data12c = cat(3,data1c,data2c);

        neworder = [8,1,7,2,6,3,5,4];
        data12c = data12c(:,:,neworder,:);
    else
        data12c = [];
        off = linspace(1,2049,5);
        data12c(:,:,1,:) = data1c(:,off(1):off(2)-1,:);
        data12c(:,:,3,:) = data1c(:,off(2):off(3)-1,:);
        data12c(:,:,5,:) = data1c(:,off(3):off(4)-1,:);
        data12c(:,:,7,:) = data1c(:,off(4):off(5)-1,:);
        data12c(:,:,2,:) = data2c(:,off(1):off(2)-1,:);
        data12c(:,:,4,:) = data2c(:,off(2):off(3)-1,:);
        data12c(:,:,6,:) = data2c(:,off(3):off(4)-1,:);
        data12c(:,:,8,:) = data2c(:,off(4):off(5)-1,:); 
    end

    if cal.correctForPSFshape == 0
        for m=(1:sys.nplanes-1)-1

            ch=sys.nplanes-m;

            img_seq_ch_fix=squeeze(data12c(:,:,ch,:));
            img_seq_ch_moving=squeeze(data12c(:,:,ch-1,:));

            [~,indcom_fix]=max(squeeze(mean(max(img_seq_ch_fix))));
            [~,indcom_mov] = max(squeeze(mean(max(img_seq_ch_moving))));

            im_fix=medfilt2(img_seq_ch_fix(:,:,indcom_fix),[3 3]);
            im_mov=medfilt2(img_seq_ch_moving(:,:,indcom_mov-2),[3 3]); % Shift from the focus plane

            % extracting the center of gravities of the beads
            out=struct;
            out.ru=5;
            sys.bg=1;
            out.nh=0;

            out=hriSegmentation(double(im_fix),cal.bgth,cal.logsize,out);
            out=hriFilterSegments(double(im_fix),cal.aupl,cal.alol,sys,out);
            cog_fix=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];
            out=hriSegmentation(double(im_mov),cal.bgth,cal.logsize,out);
            out=hriFilterSegments(double(im_mov),cal.aupl,cal.alol,sys,out);
            cog_mov=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];

            % identify corresponding beads in the two channels
            % spatial cross-correlation algorithm to determine shift of coordinates
            im_fix = im_fix(50:end-50,50:end-50);
            im_mov = im_mov(50:end-50,50:end-50);
            [mx, my] = ccrShiftEstimation(im_fix,im_mov,5);

            cog_common_fix=zeros(1,2);
            cog_common_mov=zeros(1,2);

            % identify corresponding center of gravities
            pixel_tolerance=2;
            for k=1:size(cog_fix,1)
                d=(cog_mov-repmat(cog_fix(k,:)-[my mx],[size(cog_mov,1) 1]));
                d=sqrt(d(:,1).^2+d(:,2).^2);
                [p,ind]=min(d);
                if(p<pixel_tolerance)
                    cog_common_fix(end+1,:)=cog_fix(k,:);
                    cog_common_mov(end+1,:)=cog_mov(ind,:);
                end
            end

            cog_common_fix(1,:)=[];  %first element is irrelevant (0,0)
            cog_common_mov(1,:)=[];

            d=(cog_common_mov+repmat([my mx],[size(cog_common_mov,1) 1]))-cog_common_fix;
            disp(['avg error of COG coordinates using pure displacement at pixel level: ',...
                num2str(mean(sqrt(d(:,1).^2+d(:,2).^2))),', (',num2str(mean(d(:,1))),', ',num2str(mean(d(:,2))),')']);
            cal.d{m+1}= d;
            cal.davg(m+1) = mean(sqrt(d(:,1).^2+d(:,2).^2)); % residual displacement
            cal.dmx(m+1) = mean(d(:,1));
            cal.dmy(m+1) = mean(d(:,2));
            % affine transformation
            ccm=fliplr(cog_common_mov);
            ccf=fliplr(cog_common_fix);
            tf{m+1}=cp2tform(ccm,ccf,'similarity');
            tmp = fitgeotrans(ccm,ccf,'affine');
            tf{m+1}.tdata.T = tmp.T;

            im_mov_tr=imtransform(linmap(im_mov,min(im_mov(im_mov ~= 0)),max(im_mov(:)),0,1),maketform('affine',tf{m+1}.tdata.T), 'XData', [1 size(im_fix,2)], 'YData', [1 size(im_fix,1)], ...
                'Size', size(im_fix));
            im_fix = linmap(im_fix,min(im_fix(im_fix ~= 0)),max(im_fix(:)),0,1);
            im_fix = im_fix(80:end-80,80:end-80);
            im_mov_tr = im_mov_tr(80:end-80,80:end-80);
            [mx2,my2] = ccrShiftEstimation(im_fix,im_mov_tr,5);
            tf{m+1}.tdata.T(3,1:2) = tf{m+1}.tdata.T(3,1:2) + [mx2,my2];
        end

        disp(['Integrated displacement error along z: ',num2str(sum(cal.dmx)),...
            ', ',num2str(sum(cal.dmy)),' (pixels)']);

        % ------------ correct for psf assymetry -------------
    else %  cal.correctForPSFshape == 1
        disp('Try to correct the corregistration for psf assymetry')
        tfpsf=cell(7,1);

        ch = 4;
        img_seq_ch_fix=squeeze(data12c(:,:,ch,:));
        [~,indcom_fix]=max(squeeze(mean(max(img_seq_ch_fix))));

        dx = []; dy = [];
        % sampling is 200nm but prism spacing is 350nm
        indices = indcom_fix-4:indcom_fix+3;
        indices = round((indices-indcom_fix)*350/200) + indcom_fix;
        for m=1:length(indices)-1
            ch = sys.nplanes-m;
            im_fix=medfilt2(data12c(:,:,ch+1,indcom_fix ),[3 3]);
            im_fix = im_fix(50:end-50,50:end-50);

            im_mov=medfilt2(data12c(:,:,ch, indcom_fix),[3 3]);
            im_mov = im_mov(50:end-50,50:end-50);
            % extracting the center of gravities of the beads
            out=struct;
            out.ru=5;
            sys.bg=1;
            out.nh=0;
            out=hriSegmentation(double(im_fix),cal.bgth,cal.logsize,out);
            out=hriFilterSegments(double(im_fix),cal.aupl,cal.alol,sys,out);
            cog_fix=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];
            out=hriSegmentation(double(im_mov),cal.bgth,cal.logsize,out);
            out=hriFilterSegments(double(im_mov),cal.aupl,cal.alol,sys,out);
            cog_mov=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];

            figure(1);imagesc(im_fix);hold on;plot(cog_fix(:,2),cog_fix(:,1),'rx');hold off
            figure(2);imagesc(im_mov);hold on;plot(cog_mov(:,2),cog_mov(:,1),'rx');hold off
            pause
            [mx, my] = ccrShiftEstimation(im_fix,im_mov,5);

            dx(end+1) = mx; dy(end+1) = my;
            cog_common_fix=zeros(1,2);
            cog_common_mov=zeros(1,2);

            %identify corresponding center of gravities
            pixel_tolerance=2;
            for k=1:size(cog_fix,1)
                d=(cog_mov-repmat(cog_fix(k,:)-[my mx],[size(cog_mov,1) 1]));
                d=sqrt(d(:,1).^2+d(:,2).^2);
                [p,ind]=min(d);
                if(p<pixel_tolerance)
                    cog_common_fix(end+1,:)=cog_fix(k,:);
                    cog_common_mov(end+1,:)=cog_mov(ind,:);
                end
            end

            cog_common_fix(1,:)=[];  %first element is irrelevant (0,0)
            cog_common_mov(1,:)=[];

            d=(cog_common_mov+repmat([my mx],[size(cog_common_mov,1) 1]))-cog_common_fix;
            disp(['avg error of COG coordinates using pure displacement at pixel level: ',...
                num2str(mean(sqrt(d(:,1).^2+d(:,2).^2))),', (',num2str(mean(d(:,1))),', ',num2str(mean(d(:,2))),')']);
            cal.d{m}= d;
            cal.davg(m) = mean(sqrt(d(:,1).^2+d(:,2).^2)); % residual displacement
            cal.dmx(m) = mean(d(:,1));
            cal.dmy(m) = mean(d(:,2));

            % affine transformation
            ccm=fliplr(cog_common_mov);
            ccf=fliplr(cog_common_fix);

            tmp = fitgeotrans(ccm,ccf,'affine');
            tfpsf{m}.tdata.T  = tmp.T;
        end

        for n = 1:numel(tf)
            tf{n}.tdata.T = tfpsf{n}.tdata.T;
        end

    end
    cal.tf = tf;
    %  check coregistration by applying it on beads
    coregBeads = [];
    img_seq_ch_fix=squeeze(data12c(:,:,4,:));
    [~,indcom_fix]=max(squeeze(mean(max(img_seq_ch_fix))));
    im_fix=medfilt2(squeeze(data12c(:,:,8,indcom_fix)),[3 3]);

    for m=(1:sys.nplanes-1)-1
        ch=sys.nplanes-m;
        img_seq_ch_moving=squeeze(data12c(:,:,ch-1,:));

        im_mov=medfilt2(img_seq_ch_moving(:,:,indcom_fix),[3 3]);
        T = tf{sys.nplanes-8+1}.tdata.T;
        for k = 7:-1:ch
            temp = tf{sys.nplanes-k+1}.tdata.T;
            T = T*temp;
        end

        im_mov_tr=imwarp(linmap(im_mov,min(im_mov(im_mov ~= 0)),max(im_mov(:)),0,1),imref2d(size(im_fix)),affine2d(T),...
            'OutputView',imref2d(size(im_fix)));

        coregBeads(:,:,ch-1) = linmap(im_mov_tr,0,1);
        if ch == 8
            coregBeads(:,:,8) = linmap(im_fix,min(im_fix(im_fix ~= 0)),max(im_fix(:)),0,1);
        end
    end

    % compute cropping properties from beads coregistration
    cal.mask = sum(coregBeads== 0,3) == 0;
    cal.dx = sum(cal.mask(end/2,1:end/2) == 0) + 5;
    cal.dy = sum(cal.mask(1:end/2,end/2) == 0) + 5;
    cal.Lx = sum(cal.mask(end/2,:) == 1) - 10;
    cal.Ly = sum(cal.mask(:,end/2) == 1) - 10;

    cmap = imresize(jet,[7 3],'bilinear');
    figure(1231);
    plot(cal.dmx(1),cal.dmy(1),'x','color',cmap(1,:),'linewidth',2); hold on
    for k = 2:7
        plot(cal.dmx(k),cal.dmy(k),'x','color',cmap(k,:),'linewidth',2)
    end
    for k = 1:7
        plot(cal.d{k}(:,1),cal.d{k}(:,2),'.','color',cmap(k,:));
    end
    hold off
    xlim([-1 1]); ylim([-1 1]); legend(num2str(cal.davg',3));
    title(['Integrated displacement error along z: ',num2str(sum(cal.dmx),3),...
        ', ',num2str(sum(cal.dmy),3),' (pixels)'])
    pause(0.01) % give time to Matlab to display the figure
% eof