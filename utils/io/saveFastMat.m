function saveFastMat(filename, varargin)
% fast save of large arrays to .mat files, no compression
%
% Matlab's 'save' command can be very slow when saving large arrays,
% because by default Matlab attempts to use compression. This function
% is a faster alternative with no compression.

% Append .mat if necessary
[filepath, filebase, ext] = fileparts(filename);
if isempty(ext)
    filename = fullfile(filepath, [filebase '.mat']);
end

for ii = 1:numel(varargin)
    vname = ['/' inputname(ii+1)];
    h5create(filename, vname, size(varargin{ii}), 'DataType', class(varargin{ii}));
    h5write(filename, vname, varargin{ii});
end
% eof