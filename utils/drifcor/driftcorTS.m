function stack = driftcorTS(stack,settings, frames)
% drift correction in lateral direction based on fiducial markers 
% using drift corr data from ThunderSTORM
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings
% output:   stack .. drift corrected stack

    
    data = loadjson([settings.io.imagePath,filesep,settings.io.imageName,'_drift.json']);
    
    polyy = [];
    polyx = [];

    for ii = 1:data.yFunction.n
        disp(ii)
        polyy = [polyy,polyval(data.yFunction.polynomials{ii}.coefficients,...
            ones(1,data.yFunction.knots(ii+1)-data.yFunction.knots(ii)))];

        polyx = [polyx,polyval(data.xFunction.polynomials{ii}.coefficients,...
            ones(1,data.xFunction.knots(ii+1)-data.xFunction.knots(ii)))];

    end

    polyy = [polyy,polyval(data.yFunction.polynomials{end}.coefficients,1)];
    polyx = [polyx,polyval(data.xFunction.polynomials{end}.coefficients,1)]; 
    
    drift = [polyx', polyy'];
    
    pxsize = settings.sys.pxy;
    
    [sy,sx,~] = size(stack);
    [xx, yy]=meshgrid(1:sx,1:sy);
    
    if size(stack,3) > size(drift,1)
        kmax = size(drift,1);
        stack = stack(:,:,1:kmax);
    else
        kmax = size(stack,3);
    end
    
    for k=1:kmax
        disp(['dcor: ',num2str(k)]);
        img = stack(:,:,k);

        shiftxy = drift(k,1:2);
        shiftxy = +1*shiftxy./(pxsize);

        imReg=interp2(xx,yy,double(img),xx+shiftxy(1),yy+shiftxy(2),'linear');
        imReg(isnan(imReg)) = 0;
        stack(:,:,k)=uint16(imReg);
    end

end
% eof