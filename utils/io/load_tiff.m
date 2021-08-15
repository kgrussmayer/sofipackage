function [stack,Nframes,sx,sy] = load_tiff(filename,Nframes,roi)

fileinfo = imfinfo(filename);
sx = single(fileinfo(1).Width);
sy = single(fileinfo(1).Height);

try
    frames = cellfun(@str2double,regexp(fileinfo(1).ImageDescription,'images=(\d*)','tokens'));
catch
    frames = [];
end
if isempty(frames), frames = length(fileinfo); end

if nargin < 2 || ~exist('Nframes','var') || ~any(Nframes)
    Nframes = frames;
elseif Nframes > frames
    Nframes = frames;
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
stack = zeros([sy,sx,Nframes],precision);

% read the file using low-level file I/O
fig = statusbar('Loading data...');
fid = fopen(filename,'r',machinefmt);
for n = 1:Nframes
    fig = statusbar(n/Nframes,fig);
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