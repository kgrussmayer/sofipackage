function [stack,Nframes,sy,sx] = load_bigtif(fname, Nframes,roi)

fileinfo = imfinfo(fname);
sx=fileinfo(1).Width;
sy=fileinfo(1).Height;

numFramesStr = regexp(fileinfo.ImageDescription, 'images=(\d*)', 'tokens');
frames = str2double(numFramesStr{1}{1});

if nargin < 2 || ~any(Nframes)|| isempty(frames)
    Nframes=frames;
end

if Nframes > frames
    Nframes=frames;
    disp(['Tiff file contains only ',num2str(frames),' images'])
end

if nargin <3 
roi = [];
end

% Use low-level File I/O to read the file
fp = fopen(fname , 'rb');
% The StripOffsets field provides the offset to the first strip. Based on
% the INFO for this file, each image consists of 1 strip.
fseek(fp, fileinfo.StripOffsets, 'bof');
% Assume that the image is 16-bit per pixel and is stored in big-endian format.
% Also assume that the images are stored one after the other.
  
stack=zeros(sy,sx,Nframes,'uint16');
fig=statusbar('Loading data...'); 
for ii = 1:Nframes
    fig=statusbar(ii/Nframes,fig);
    stack(:,:,ii) = fread(fp, [sx sy], 'uint16', 0, 'ieee-be')';
end
fclose(fp);
delete(fig);

if ~isempty(roi)
stack = stack(roi(1,:),roi(2,:));
end


