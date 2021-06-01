function stack = driftcor(stack,settings)
% drift correction in lateral direction based on fiducial markers 
% requires additonal mat file with x y shifts in each time frame
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings
% output:   stack .. drift corrected stack

    drift = load([settings.io.imageFile,settings.dcor.tag,'.mat']);
    drift = drift.drift;
%     drift = drift_corr;
%      drift = drift;
%     pxsize = 104.8; % projected pixel size in nm (Hendrik's setup)
%     pxsize = 16000/166.6; % projecte pixel size in nm (sofi setup)
    
    pxsize = settings.sys.pxy;
    
    [sy,sx,~] = size(stack);
    [xx, yy]=meshgrid(1:sx,1:sy);
    
    if size(stack,3)-1 <size(drift,1)
        kmax = size(stack,3)-1;
    else
        kmax = size(drift,1);
    end
    
    for k=1:kmax
        disp(['dcor: ',num2str(k)]);
        img = stack(:,:,k+1);

        shiftxy = drift(k,2:3);
        shiftxy = +shiftxy./(pxsize);

        imReg=interp2(xx,yy,double(img),xx+shiftxy(1),yy+shiftxy(2),'linear');
        imReg(isnan(imReg)) = 0;
        stack(:,:,k+1)=uint16(imReg);
    end

end