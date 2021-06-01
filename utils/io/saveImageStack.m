function saveImageStack(stack,outputPath,fileName,tag,bits)

if nargin <4
    tag = '';
    bits = 8;
end

if nargin <5
    bits = 8;
end

nowID = datestr(now,'ddmmyyyyHHMMSS');

% outputPath = [outputPath, filesep,fileName];

if (~exist(outputPath,'dir'))
    mkdir(outputPath);
end

% outputFile = [outputPath,filesep,fileName,'_',tag,nowID]; % create unique name with date
outputFile = [outputPath,filesep,fileName,tag];

stack = stack./max(stack(:));

if bits == 16; stack = uint16(stack.*2^16);end;

fig = statusbar('Saving image stack...');
frames = size(stack,3);
for jj = 1:frames
    fig = statusbar(jj/frames,fig);
    imwrite(stack(:,:,jj),[outputFile,'.tif'],'WriteMode', 'append',  'Compression','none');
    pause(0.1);
end
delete(fig);

disp(['figure was succesfully saved into: ',outputFile]);

end