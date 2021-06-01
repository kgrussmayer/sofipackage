function saveFigure(fh,outputPath,fileName,format,addID)
% saves the figure fh according to outputPath and fileName
% format .. type of the output file
% i.e. pdf, png, tif, fig etc. 

% Example of use:
% % save the current figure in tif file
% fh = gcf;saveFigure(fh,outputPath,'nameOfTheFigure','tif'); 

if nargin <4
    format='fig';
    addID = 0;
end

if nargin <5
    addID = 0;
end

% outputPath = [outputPath, filesep,fileName];

% create the output directory if it does not exists
if (~exist(outputPath,'dir'))
    mkdir(outputPath);
end

if addID == 1
    nowID = datestr(now,'mmmm-dd-yyyy_HHMMSS');
    % add a unique identifier to the file name
    outputFile = [outputPath,filesep,fileName,'_',nowID];
else
    outputFile = [outputPath,filesep,fileName];
end

set(fh, 'PaperPositionMode','auto');
saveas(fh,[outputFile,'.',format]);
%disp(['figure was succesfully saved into: ',outputFile]);

end