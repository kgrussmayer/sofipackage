%
% Main for higher order SOFI using biplane data
%

% Load all the processing settings from config file
% (see experiments)

addpath(genpath('utils'));
settings.method = 'sofibp';

%% loading data

data1 = load_tifFile([settings.io.imagePath,filesep, ...
    settings.io.imageName, settings.io.fext],settings.sys.sub);
data2 = load_tifFile([settings.io.imagePath2,filesep, ...
    settings.io.imageName2, settings.io.fext],settings.sys.sub);

data2 = flip(data2, 2);

[~,~,frames1] = size(data1);
[~,~,frames2] = size(data2);
frames = min([frames1,frames2]);

% apply channel weights
ch_weights = settings.sys.ch_weights./max(settings.sys.ch_weights (:));
ch_weights = 1./ch_weights;

data1 = data1 * ch_weights(1);
data2 = data2 * ch_weights(2);
%% run calibration  
disp('calibration...');

if settings.cal.beads == 1
    cnames = getnames(settings.io.pnc,'*calib.tmat');
    if ~isempty(cnames)
        assert(length(cnames)==1, 'Multiple calibration files detected.')
        cal = load(cnames{1});
    else
        cal = runcalib_biplane(settings.sys,settings.cal,settings.io); 
    end
end
    
%% register planes and calculate 3D cumulants 
sofizt = cell(1,numel(settings.sys.orders));
if isempty(settings.sys.start); settings.sys.start = 1; end

if settings.sys.subseqlength > frames
    settings.sys.subseqlength = frames;
    disp('Subsequence length is larger than available frames.')
    disp(['max subseqlength ', int2str(frames), ' will be used'])
end

for jj = settings.sys.start:settings.sys.subseqlength:...
        ceil((frames-settings.sys.start+1)/settings.sys.subseqlength)*settings.sys.subseqlength
    disp(jj)
    if (jj+settings.sys.subseqlength-1) > size(data1,3) && (size(data1,3)-jj > 500)
        st = min([size(data1,3),min(size(data2,3))]);
        data1sub = data1(:,:,jj:st);
        data2sub = data2(:,:,jj:st);
    else
        data1sub = data1(:,:,jj:jj+settings.sys.subseqlength-1);
        data2sub = data2(:,:,jj:jj+settings.sys.subseqlength-1);
    end
    sofiz = cell(1,numel(settings.sys.orders));

    if settings.cal.beads == 1
        disp('warping...');
        data1subw = imwarp(data1sub, cal.tf{1},'OutputView',imref2d(size(data2sub)));
    else
        data1subw = data1sub;
    end

    % refine registration
    disp('refine registration...');
    [sofi1]=sofiCumulants2D(data1subw,[],[],[],2);
    [sofi2]=sofiCumulants2D(data2sub,[],[],[],2);
    sofi1 = sofi1{2};
    sofi2 = sofi2{2};
    
    [mx, my] = ccrShiftEstimation(sofi1,sofi2,2); 
    shiftxy = [mx, my];

    [sy,sx] = size(sofi1);
    [xx, yy]=meshgrid(1:sx,1:sy);
    sofi1w=interp2(xx,yy,double(sofi1),xx+shiftxy(1),yy+shiftxy(2),'linear');
    sofi1w(isnan(sofi1w)) = 0;
 
    shiftxy = [mx, my]./2;
    settings.sys.shiftxy = shiftxy;
    [sy,sx,sz] = size(data1subw);
    [xx, yy]=meshgrid(1:sx,1:sy);

    for k=1:size(data1subw,3)
        disp(['dcor: ',num2str(k)]);
        img = data1subw(:,:,k);

        imReg=interp2(xx,yy,double(img),xx+shiftxy(1),yy+shiftxy(2),'linear');
        imReg(isnan(imReg)) = 0;
        data1subw(:,:,k)=imReg;
        
    end

    % calculate SOFI cumulants

        % 3D SOFI
        [sofi3D]=sofiCumulants3D(permute...
            (cat(4,data1subw,data2sub),[1 2 4 3]),[],[],[],settings.sys.orders);

        % 2D SOFI for the previously transformed plane
        [sofi2D]=sofiCumulants2D(data1sub,[],[],[],settings.sys.orders);

    % reassemble the 3D stack
    for order =  settings.sys.orders
        sofiz{order} = cat(3,sofi2D{order},sofi3D{order}(:,:,2:end));

        sofizt{order} = cat(4,sofizt{order},sofiz{order});   
    end
end
   
%% final alignment of SOFI stack
% realign the whole 3D stack
% warp on 1 SOFI plane, shifts in between all planes in the stack
disp('final realignment...');

for ii = settings.sys.orders

    sofi = sofizt{ii};
    [sy,sx,sz,st] = size(sofi);
    [xx, yy]=meshgrid(1:sx,1:sy);
    shiftxy = settings.sys.shiftxy*ii;

    for ff = 1:st
        if settings.cal.beads == 1
            img = imwarp(sofi(:,:,1,ff), cal.tf{ii},'OutputView',imref2d(size(sofi(:,:,1,ff))));
            sofi(:,:,1,ff)=interp2(xx,yy,img,xx+shiftxy(1),yy+shiftxy(2),'linear');
            sofi(isnan(sofi))=0; 
        else
            sofi(:,:,1,ff)=interp2(xx,yy,sofi(:,:,1,ff),xx+shiftxy(1),yy+shiftxy(2),'linear');
            sofi(isnan(sofi))=0; 
        end
    end
    sofizt{ii} = sofi; 

end
      
%% Perform flattening 
% distance factor correction in all spatial dimensions
disp('Flattening');
sofizt=sofiAllFlatten3D(sofizt,settings.sys.orders);

%% Deconvolution
disp('deconvolution and linearization');
dec = settings.dec;
dec.avgn = 1;

if dec.avgn > 1
sofizt=cellfun(@(x)squeeze(mean(reshape(x,size(x,1),size(x,2),size(x,3),dec.avgn,[]),4)),sofizt,'UniformOutput',0);
end

% apply median filtering if it is set on
if dec.medfilt == 1
    % 2D median filtering
    for n=1:size(sofizt{2},3)
        for m=1:size(sofizt{2},4)
            sofizt{2}(:,:,n,m)=medfilt2(sofizt{2}(:,:,n,m),[3 3]);
        end
    end
end

[sofid]=sofiLinearize3Dpad(sofizt,dec);

%% Save results

% create the output folder
if (~exist(settings.io.outputpath,'dir'))
    mkdir(settings.io.outputpath);
end

% export output stacks as tif files
for ii = settings.sys.orders %skip firrst order
    
    sofi_l = mean(sofid{ii},4);
    saveImageStack(sofi_l,settings.io.outputpath,...
        [settings.io.imageName,'_ord_',num2str(ii),'_lucy'],[],16);
    
    sofi_c = mean(sofizt{ii},4);
    saveImageStack(sofi_c,settings.io.outputpath,...
        [settings.io.imageName,'_ord_',num2str(ii),'_sofic16'],[],16);
    
end   

