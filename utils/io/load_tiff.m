function [stack,frames,sx,sy] = load_tiff(filename,n_frames,roi)
% Load a stack of images from a multiple tif file (supports both standard
% and big tiff files).
% file_name contains path to the file to be loaded. Nframes is a number of
% frames which are going to be loaded from the multiple tif file. Optinally
% only roi (region of interest) can be loaded.
%
% Inputs:
% file_name     path to the multiple tif file
% n_frames      number of frames to be loaded
% roi           region of interest to be loaded
%
% Outputs:
% stack         stack of images (rows, columns, frames)
% frames        number of frames in the loaded stack
% sx            image width (number of columns)
% sy            image height (number of rows)

fileinfo = imfinfo(filename);
sx = single(fileinfo(1).Width);
sy = single(fileinfo(1).Height);

try
    frames = cellfun(@str2double,regexp(fileinfo(1).ImageDescription,'images=(\d*)','tokens'));
catch
    frames = [];
end
if isempty(frames), frames = length(fileinfo); end

if nargin < 2 || ~exist('n_frames','var') || ~any(n_frames)
    n_frames = frames;
elseif n_frames > frames
    n_frames = frames;
    disp(['Tiff file contains only ',num2str(frames),' images'])
end

if nargin < 3
    roi = [];
end

switch fileinfo(1).ByteOrder
    case 'little-endian'
        machinefmt = 'l';
    case 'big-endian'
        machinefmt = 'b';
    otherwise
        machinefmt = 'n';   % native system byte ordering
end

switch fileinfo(1).BitDepth
    case 8
        precision = 'uint8=>uint16';
    case 16
        precision = 'uint16';
    case 32
        precision = 'uint32';
    otherwise
        precision = 'double';
end
stack = zeros([sy,sx,n_frames],precision);

% read the file using low-level file I/O
fig = statusbar('Loading data...');
fid = fopen(filename,'r',machinefmt);
for n = 1:n_frames
    fig = statusbar(n/n_frames,fig);
    fseek(fid,fileinfo(n).StripOffsets(1),'bof');
    stack(:,:,n) = fread(fid,[sx,sy],precision,0,machinefmt)';
end
fclose(fid);
delete(fig);

if ~isempty(roi)
    stack = stack(roi(1,:),roi(2,:),:);
end
[sy,sx,~] = size(stack);
% eof