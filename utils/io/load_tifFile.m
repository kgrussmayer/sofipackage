function stack = load_tifFile(fname, Nframes,roi)
% loads a stack of images from a multiple tif file,which is specified by its path and name
% fname is a string containing path and name of the file to be loaded
% Nframes is a number of frames which are going to be loaded from the multiple tif file

fileinfo = imfinfo(fname);

if nargin < 3
    roi = [];
end

if nargin < 2 || ~any(Nframes) || isempty(Nframes)
    Nframes=length(fileinfo);
end

if Nframes > length(fileinfo)
        Nframes=length(fileinfo);
        disp(['Tiff file contains only ',num2str(Nframes),' images'])
end

% stack =  zeros(fileinfo(1).Height,fileinfo(1).Width,Nframes); %prelocate memory 
fig=statusbar('Loading data...');

for ii=1:Nframes
%     disp(['Loading.. ',num2str(ii)]);
    fig=statusbar(ii/Nframes,fig);
    temp=uint16(imread(fname,ii,'info',fileinfo));
    if size(temp,3) > 1 % RGB image
        temp = rgb2gray(temp);
    end
    
    if any(roi)
       stack(:,:,ii) = temp(roi(1,:),roi(2,:));
    else
       stack(:,:,ii) = temp; 
    end
end
delete(fig);