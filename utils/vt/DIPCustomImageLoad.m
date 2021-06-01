function [ Stack,StackO] = DIPCustomImageLoad( S )
% S ... Structure With Input Parameters
% Stack ... Requested Image
% StackO ... Original Image

% EXTRACT INPUT PARAMS 
ImagePath = S.ImagePath;
ImageNames = S.ImageNames;
ImName = S.ImName;
ImStart = S.ImStart;
NumImages = S.NumImages;
PixLims = S.PixLims;
Nmax = S.NMax;

    
% LOAD REAL IMAGE DATA
disp('LOADING IMAGE ...');
ImageInfo = imfinfo(horzcat(ImagePath,ImName,'.tif'));
Slice = ImageInfo(1);
N = length(ImageInfo);
if(N < NumImages + ImStart)
    error(horzcat('NOT ENOUGH IMAGE FRAMES (Requested / Available): ',num2str(ImStart + NumImages),' / ',num2str(N)));
end

% GET PREALLOCATING DIMS
MaxWidth = Slice.Width;
MaxHeight = Slice.Height;
if(MaxWidth > PixLims || MaxHeight > PixLims)
    XStart = floor(MaxWidth/2 - PixLims/2);
    XStop = floor(MaxWidth/2 + PixLims/2)-1;
    YStart = floor(MaxHeight/2 - PixLims/2);
    YStop = floor(MaxHeight/2 + PixLims/2)-1;
    StackO = single(zeros(PixLims,PixLims,min(N,Nmax)));
else
    StackO = single(zeros(MaxHeight,MaxWidth,min(N,Nmax)));
end

% LOAD IMAGE
for i = 1:min(N,Nmax)
    Slice = im2double(imread(horzcat(ImagePath,ImName,'.tif'),i,'info',ImageInfo));
    if(size(Slice,3) > 1)
        Slice = rgb2gray(Slice);
    end

    if(MaxWidth > PixLims || MaxHeight > PixLims)
        Slice = Slice(YStart:YStop,XStart:XStop);
    end
    StackO(:,:,i) = Slice;
end

Stack = StackO(:,:,ImStart:ImStart + NumImages-1);
    


end








