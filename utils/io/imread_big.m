%Tristan Ursell
%Read large image stack (TIFF)
%May 2019
%
% This function can load image stacks larger than 4GB which is the typical
% limitation on image files.  Designed to work with single-channel 
% uncompressed TIFF stacks.  Also works for files smaller than 4GB.
%
% [stack_out,Nframes]= imread_big(stack_name);
% [stack_out,Nframes]= imread_big(stack_name,[i j]);
%
% stack_name = the path and file name to the image stack
%
% [i j] = optional frame number range to load (i = j is allowed)
%
% Nframes = number of frames as determined by total file size divided by
% estimated size of each frame data block
%
% stack_out = the output image stack of size [M N Nframes], where M x N is
% the size of each image.
%

function [stack_out,Nframes] = imread_big(stack_name,varargin)

%get data block size
info1 = imfinfo(stack_name);
stripOffset = info1(1).StripOffsets;
stripByteCounts = info1(1).StripByteCounts;

%get image size
sz_x=info1(1).Width;
sz_y=info1(1).Height;
if length(info1)<2
    Nframes=floor(info1(1).FileSize/stripByteCounts);
else
    Nframes=length(info1);
end

%read in arguments
if nargin>1
    temp1=varargin{1};
    ilow=temp1(1);
    ihigh=temp1(2);
else
    ilow=1;
    ihigh=Nframes;
end

%load file id tag
fID = fopen (stack_name, 'r');

if info1(1).BitDepth==32
    stack_out = zeros([sz_y sz_x ihigh-ilow+1],'uint32');
elseif info1(1).BitDepth==16
    stack_out = zeros([sz_y sz_x ihigh-ilow+1],'uint16');
else
    stack_out = zeros([sz_y sz_x ihigh-ilow+1],'uint8');
end

start_point = stripOffset(1) + (0:1:(Nframes-1)).*stripByteCounts + 1;
%start_point = stripOffset(1) + (0:1:(Nframes-1))*stripByteCounts + 1;

q=0;
for i = ilow:ihigh
    %fprintf ('loading image ... %d\n', i);
    fseek (fID, start_point(i), 'bof');
    
    if info1(1).BitDepth==32
        A = fread(fID, [sz_x sz_y], 'uint32=>uint32');
    elseif info1(1).BitDepth==16
        A = fread(fID, [sz_x sz_y], 'uint16=>uint16');
    else
        A = fread(fID, [sz_x sz_y], 'uint8=>uint8');
    end
    
    q=q+1;
    try
        stack_out(:,:,q) = A';
    catch
        disp(['Terminated at frame ' num2str(i) '.'])
        break
    end
end

%stack_out=stack_out(:,:,1:q);

Nframes=q;
fclose(fID);




