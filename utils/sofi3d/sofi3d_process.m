function [sofizt,sys,cal] = SOFI3D_process(sys,cal,io,pn1) 
% Perform all steps of the multiplane SOFI algorithm
% author: Tomas Lukes, tomas.lukes@epfl.ch, parts of the code taken from Stefan's codes
% MATLAB version 2012a

% note: This script calculates cumulants for 3D images acquired by 2 digital cameras
% 1) 2D raw cumulants are calculated for each z-plane
% 2) Calibration files (images of fluorescent beads) are loaded and transformation matrices between consecutive 2D raw cumulants are calculated 
% 3) 2D raw cumulants are registered and 3D cumulants are calculated

%% Initiate
if strcmp(sys.fext,'.bin')
    sys.totlength = loadbinfileLength([pn1 filesep],io.fn1); %[pn1 filesep]
elseif strcmp(sys.fext,'.dat')
    temp = dir(pn1);
    sys.totlength = temp.bytes/(2*2*2048*600);
else
    info = imfinfo(pn1);
    h = info.Height; w = info.Width; fsize = info.FileSize;
    sys.totlength = floor(fsize./(2*h.*w)/2);
end
temp.path = pn1;
if isempty(sys.start), sys.start = 1; end
if isempty(sys.nss) 
    % take into account overlap of subsequences if sys.subseqstep is smaller
    % than sys.subseqlength
    Nss=floor((sys.totlength-sys.start+1-sys.subseqlength)/sys.subseqstep+1);
else
    Nss = sys.nss;
end

sys.size(1) = 512;
sys.size(2) = 512;
sys.camorigin = [1,1,1,1,2,2,2,2];

%% Calibration - calculate transformation matrices between 2D raw cumulants in different z-planes                   

if 1%~isfield(cal,'tf')
    disp('Calibration calculation')
    cal = runcalib(sys,cal,io.pnc); 
    if cal.saveCal == 1
       saveas(1231,[io.pn3,filesep,'coregistration error'],'png');
       h = guidata(1230); % grab coregistered beads data
       writeTIFF(h.data,[io.pnc,filesep,getID(5),'_coregistered_beads']);
       close(1230);
       pause(0.1);
       disp('Coregistration results saved in calib folder')
    end
end


%% 3D cumulants calculation

  disp(['Data : ',pn1])
    
    if strcmp(sys.fext,'.bin')
    	disp([pn1 filesep 'data1.bin']);
        fid1=fopen([pn1 filesep 'data1.bin']);
        fid2=fopen([pn1 filesep 'data2.bin']);
    elseif strcmp(sys.fext,'.dat')
        fid = fopen(pn1);
    else
    	disp('loading data')
        dataTemp = imread_big(pn1);
        % flip camera 2
        for k = 2:2:size(dataTemp,3)
        	dataTemp(:,:,k) = fliplr(dataTemp(:,:,k));
        end
    end
    
    sofizt = cell(1,numel(sys.orders));

    for ns=1:Nss
        disp(['Processing subsequence ',num2str(ns),' of ',num2str(Nss)])
        sofiz = cell(1,numel(sys.orders));
        
        if strcmp(sys.fext,'.bin')
            tic
        % allow overlap of subsequences 
            fr=sys.start+(ns-1)*sys.subseqstep-1+(1:sys.subseqlength);
        
        % move to starting frame in data1
            fseek(fid1, 2*(fr(1)-1).*sys.size(1).*sys.size(2).*round(sys.nplanes/max(sys.camorigin)), 'bof');
        % move to starting frame in data2
            fseek(fid2, 2*(fr(1)-1).*sys.size(1).*sys.size(2).*round(sys.nplanes/max(sys.camorigin)), 'bof');
         
            data1=uint16(zeros(sys.size(1),sys.size(2).*round(sys.nplanes/max(sys.camorigin)),numel(fr)));
            data2=data1;
        
            for nf=fr
                data1(:,:,-fr(1)+1+nf)=fread(fid1,[sys.size(1),sys.size(2).*round(sys.nplanes/max(sys.camorigin))],'*uint16');
                data2(:,:,-fr(1)+1+nf)=fread(fid2,[sys.size(1),sys.size(2).*round(sys.nplanes/max(sys.camorigin))],'*uint16');
            end
        elseif strcmp(sys.fext,'.dat')
            % read data from Labview software
            disp('Loading Labview data')
            t0 = tic;
            fr=2*(sys.start+(ns-1)*sys.subseqstep)-1+(1:2*sys.subseqlength);
            fseek(fid, 2*(fr(1)-1)*2048*600,'bof');
            data1 = uint16(zeros(600,2048,sys.subseqlength));
            data2 = data1;
            count = 1; id = 1;
            for nf = fr
                if mod(count,2) == 1
                    temp = fread(fid,[2048,600],'uint16','ieee-be');
                    data1(:,:,id) = temp';
                else
                    temp = fliplr(fread(fid,[2048,600],'uint16','ieee-be'));
                    data2(:,:,id) = rot90(temp',2);
                    id = id+1;
                end
                count = count +1;
            end
            disp(['Loading time : ',num2str(toc(t0)),' [s]'])
            % TODO: read data size from .info file
        else
            % get subsequences from full stack
            fr=sys.start+(ns-1)*sys.subseqstep-1+(1:sys.subseqlength);
            data1 = dataTemp(:,:,2*fr(1):2:2*fr(1)+2*sys.subseqlength-1);
            data2 = dataTemp(:,:,2*fr(1)+1:2:2*fr(1)+2*sys.subseqlength);
            % TODO: modify imread_big to load only selected frames
        end

        %%% Reorder planes for data using new prism
        if strcmp(sys.fext,'.bin')
            data1 = reorderData(data1,1,0,sys);
            data2 = reorderData(data2,2,1,sys);
        
            [sy,sx,st] = size(data1);
            data1 = reshape(data1,sy,sx/4,4,st);
            data2 = reshape(data2,sy,sx/4,4,st);
            data12 = cat(3,data1,data2);
            % TODO order is setup specific - selection based on imaging
            % system (has to be the same as the one used in runcalib
%             neworder = [2,7,1,8,4,5,3,6]; % plane order of labview setup
%             (dat files)
            neworder = [8,1,7,2,6,3,5,4]; % plane order of sofi new prism
%             setup (bin files)
            data12 = data12(:,:,neworder,:);
        
%             data1 = data12(:,:,1:4,:);
%             data1 = reshape(data1,sy,sx,st);
%             data2 = data12(:,:,5:8,:);
%             data2 = reshape(data2,sy,sx,st);

        else
        % crop data & reorder planes for Labview and Micromanager data
            disp('Reordering planes')
            data12 = zeros(size(data1,1),512,8,sys.subseqlength,'uint16');
            off = linspace(1,2049,5);
            data12(:,:,1,:) = data1(:,off(1):off(2)-1,:);
            data12(:,:,3,:) = data1(:,off(2):off(3)-1,:);
            data12(:,:,5,:) = data1(:,off(3):off(4)-1,:);
            data12(:,:,7,:) = data1(:,off(4):off(5)-1,:);
            clear data1;
            data12(:,:,2,:) = data2(:,off(1):off(2)-1,:);
            data12(:,:,4,:) = data2(:,off(2):off(3)-1,:);
            data12(:,:,6,:) = data2(:,off(3):off(4)-1,:);
            data12(:,:,8,:) = data2(:,off(4):off(5)-1,:);
            clear data2;
        end

        disp('Computing 3D cumulants')
        for ch=sys.nplanes:-1:2
            fprintf('%d, ',ch);

            T=cal.tf{sys.nplanes-ch+1}.tdata.T;

            % imwarp is faster than imtransform
            im_mov_tr=imwarp(squeeze(data12(:,:,ch-1,:)),imref2d(size(squeeze(data12(:,:,ch-1,:)))),affine2d(T),...
                'OutputView',imref2d(size(squeeze(data12(:,:,ch-1,:)))));
            if ns == 1
                t1 = linmap(std(single(squeeze(data12(:,:,ch  ,:))),[],3),0,1);
                t2 = linmap(std(single(squeeze(data12(:,:,ch-1,:))),[],3),0,1);
                t2 = imwarp(t2,imref2d(size(t2)),affine2d(T),...
                    'OutputView',imref2d(size(t2)));
                t(:,:,:,sys.nplanes-ch+1) = mergeToRgb(t1,t2,t1./2); % magenta/green overlay
                figure(50);imagesc(t(:,:,:,end))
            end
            [sofi]=sofiCumulants3D(permute(cat(4,squeeze(data12(:,:,ch,:)),im_mov_tr),[1 2 4 3]),[],[],[],sys.orders);

            for order=sys.orders
                sofiz{order} = cat(3,sofiz{order},sofi{order}(:,:,1:end-1));
            end
            
        end
        if ns == 1
            writeRGBTIFF(uint8(255.*t),[io.pnc,filesep,getID(5),'_raw_planes_coregistration_',num2str(ch)])
        end
        disp(' ')
        [sofi]=sofiCumulants2D(squeeze(data12(:,:,ch-1,:)),[],[],[],sys.orders);
        
        
        for order=sys.orders
        sofiz{order} = cat(3,sofiz{order},sofi{order});
        sofizt{order} = cat(4,sofizt{order},sofiz{order});
        end
        
    end

    clear data12;
    clear im_mov_tr;
    
    if strcmp(sys.fext,'.bin')
        fclose(fid1);
        fclose(fid2);
    end

