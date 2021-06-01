clear all
addpath(genpath('utils'));

%% select multiple matfiles

[fileName, pathName, filterIndex] = uigetfile({'*.mat';}, 'Select file(s)', 'MultiSelect', 'on');
if iscell(fileName)
    nbfiles = length(fileName);
elseif fileName ~= 0
    nbfiles = 1;
else
    nbfiles = 0;
end

%% load molecular parameters from all selected mat files

for ii = 1:length(fileName)
    matfile_data = load([pathName,filesep,fileName{ii}]);
    bsofi_all{ii} = matfile_data.bsofi;
    ratio_all{ii} = matfile_data.ratio;
    density_all{ii} = matfile_data.density;
    brightness_all{ii} = matfile_data.brightness;
end
 
%% check array size 
cells_size = cellfun(@size, bsofi_all, 'UniformOutput', false);
cells_size = reshape(cell2mat(cells_size),2,[]);

if all(cells_size(1,:) ~= cells_size(1,1))
    disp('WARNING: input arrays do not have the same number of rows')
end

if all(cells_size(2,:) ~= cells_size(2,1))
    disp('WARNING: input arrays do not have the same number of columns')
end

rows = cells_size(1,1);
columns = cells_size(2,1);
frames = length(bsofi_all);
bsofi_all = reshape(cell2mat(bsofi_all), rows, columns, []);
ratio_all = reshape(cell2mat(ratio_all), rows, columns, []);
density_all = reshape(cell2mat(density_all), rows, columns, []);
brightness_all = reshape(cell2mat(brightness_all), rows, columns, []);

%% Merge images with drif correction
im1 = bsofi_all(:,:,1);  % reference image

drifts = zeros([frames, 2]);
for ii = 1:frames
    disp(ii)
    im2 = bsofi_all(:,:,ii); 
    [drift_x, drift_y] = ccrShiftEstimation(im1, im2, 2);
     drifts(ii, :) = [drift_x, drift_y];
end

bsofi_all = interpolateDriftcor(bsofi_all, drifts);
ratio_all = interpolateDriftcor(ratio_all, drifts);
density_all = interpolateDriftcor(density_all, drifts);
brightness_all = interpolateDriftcor(brightness_all, drifts);

%% crop edges
margin_x = 35;
margin_y = 35;
bsofi_all = bsofi_all(margin_y+1:end-margin_y, margin_x+1:end-margin_x, :);
ratio_all = ratio_all(margin_y+1:end-margin_y, margin_x+1:end-margin_x, :);
density_all = density_all(margin_y+1:end-margin_y, margin_x+1:end-margin_x, :);
brightness_all = brightness_all(margin_y+1:end-margin_y, margin_x+1:end-margin_x, :);
%%
alpha = mean(bsofi_all,3);
alpha = nperc(medfilt2(alpha,[3 3]),0.2);
alpha_max = max(alpha(:));
alpha = alpha./max(alpha(:));
alpha(alpha>=0.06) = 1;
alpha(alpha<0.06) = 0;

ratio = mean(ratio_all,3);
bright = mean(brightness_all,3);
density = sum(density_all,3);
bsofi = mean(bsofi_all,3);
%%
settings.io.outputpath = pathName;
settings.io.imageName = 'MTdnaP_650_2_MMStack_'; 

figure,
imshow(mean(bsofi_all,3),[0 0.05]);
title('bSofi')
saveImageStack(bsofi,settings.io.outputpath,[settings.io.imageName,'bsofi'],[],16);

%%
figure,
imshow(brightness,[0 1000])
title('Brightness (dn)');colormap('jet');colorbar;
saveImageStack(brightness,settings.io.outputpath,[settings.io.imageName,'brightness'],[],16);
colorbar;
%%
brightness = mean(brightness_all,3);
alpha = zeros(size(brightness));
alpha(brightness>50) = 1;  % 500 for background subtraction
tmp = alpha.*mean(ratio_all,3);
tmp(tmp<=0) = [];

figure,
imshow(alpha.*ratio,[0 0.06])
colorbar;
title('On-time ratio (-)');colormap('jet');colorbar;
saveImageStack(alpha*ratio,settings.io.outputpath,[settings.io.imageName,'ratio'],[],16);

density_saturated = nperc(density,0.2);
figure,
imshow(alpha.*density_saturated,[0 10])
title('Density (emitters/px^2)');colormap('jet');colorbar;
saveImageStack(alpha*density_saturated,settings.io.outputpath,[settings.io.imageName,'density'],[],16);
