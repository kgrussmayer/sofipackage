
function writeTIFuint(data, tiffname,res)

if nargin < 3
    res = [1 1]; % save generic resolution if none provided
end

%% save to tiff while keeping floating value

% writes data as a multi-channel TIFF with single prec. float pixels
filename = [tiffname '.tif'];
maxCount = 10;
count = 1;
while exist(filename, 'file') == 2
    filename = [tiffname '_' num2str(count) '.tif'];
    count = count +1;
    if count == maxCount
        break;
    end
end
    
    for k = 1:size(data,3)
        t = Tiff(filename, 'a');
        tagstruct.ImageLength = size(data, 1);
        tagstruct.ImageWidth = size(data, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.XResolution = res(1);
        tagstruct.YResolution = res(2); % YResolution encode z-res
%         tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = 1;
        tagstruct.BitsPerSample =  16;                        % float data
        tagstruct.SamplesPerPixel = 1;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        t.setTag(tagstruct);
        t.write(uint16(data(:,:,k)));
        t.close();
        pause(0.01)
    end

end