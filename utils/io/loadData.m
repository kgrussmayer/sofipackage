function [im,path] = loadData(path)
if nargin < 1
    [fname,pname] = uigetfile('*.*');
    path = [pname,filesep,fname];
end

ind = strfind(path,'.');
tag = path(ind(end)+1:end);
if strcmp(tag,'tif') || strcmp(tag,'tiff')
    info = imfinfo(path);
    w = info(1).Width; h = info(1).Height;
    bd = info(1).BitDepth/8;
    nf = floor(info(1).FileSize/(bd*w*h));
    keepReading = 1; k = 1;
    im = zeros(h,w,nf);
    warning('off')
    while keepReading
        try
            im(:,:,k) = imread(path,k);
            k = k+1;
        catch
            keepReading = 0;
            disp('Finished reading... ')
            disp(['Stack size : ',num2str(size(im))])
        end
    end
    warning('on')
elseif strcmp(tag,'dat')
    
    warning('off')
    fin =fopen(path,'r');
    info = dir(path);
    
    pathInfo = strrep(path,'.dat','.info');
    infoLabview = getInfoLabview(pathInfo);
    row = str2double(infoLabview.ROI_CAM0.VWid);
    col = str2double(infoLabview.ROI_CAM0.HWid);
    if row == 0
        row = 600;
        col = 2048;
    end
    nframes = info.bytes/(2*row*col);
    
    I=fread(fin,row*col*nframes,'uint16=>uint16');
    im=reshape(swapbytes(I),col,row,nframes);
    im = permute(im,[2 1 3]);
    fclose(fin);
    disp('Finished reading... ')
    disp(['Stack size : ',num2str(size(im))])
    
elseif strcmp(tag,'bin')
    warning('off')
    fin =fopen(path,'r');
    info = dir(path);
    
    pathInfo = strrep(path,'.bin','.info');
    infoLabview = getInfoLabview(pathInfo);
    row = str2double(infoLabview.ROI_CAM0.VWid);
    col = str2double(infoLabview.ROI_CAM0.HWid);
    if row == 0
        row = 600;
        col = 2048;
    end
    nframes = info.bytes/(2*row*col);
    
    I=fread(fin,row*col*nframes,'uint16=>uint16');
    im=reshape(swapbytes(I),col,row,nframes);
    im = permute(im,[2 1 3]);
    fclose(fin);
    disp('Finished reading... ')
    disp(['Stack size : ',num2str(size(im))])
end
% eof