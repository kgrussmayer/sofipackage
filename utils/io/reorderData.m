function data = reorderData(data,camid,doflip,sys,doweights)
% reorder planes, correct differences in intensities
% camid .. number of the camera {1,2}
% channel average intensities - according to beads calibration
%     ch_weights = [166.6, 151.2, 139.2, 154.5;
%     orange           141.75, 147.2, 148.1, 181];
% red established Dec. 2016 average of 6 scans [0.88, 0.93, 0.73, 1; 0.88, 0.82, 0.80, 0.94]
% orange established Dec. 2016 average of 4 scans [0.93, 0.88, 0.80, 1.00; 0.84, 0.81, 0.85, 0.99]

% ch order: cam1; cam2 -> 2,4,6,8 ; 7,5,3,1

if any(reshape(sys.ch_weights==0,1,[]))
    disp('WARNING CHANNEL WEIGTS CONTAIN ZERO VALUES!')
    return;
end

if nargin < 5
    doweights = 1;
end

ch_weights = sys.ch_weights ./max(sys.ch_weights (:));
ch_weights = 1./ch_weights;

[sy,sx,st] = size(data);
imsr = reshape(data,sy,sx/4,4,st); % reshape to get y,x,z,t


for ii = 1:4
    % flip images (left to right)
    if doflip == 1
        imsr(:,:,ii,:) = flip(imsr(:,:,ii,:),2); % for every z plane, flip image along x
    end
    
    if doweights ==1
        % correct different intensity levels
        imsr(:,:,ii,:) = imsr(:,:,ii,:)*ch_weights(camid,ii);
    end
end

data = squeeze(reshape(imsr,sy,sx,1,st));
% eof