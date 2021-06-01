% Run Multicolor SOFI. Tested with Win7, Win10, Matlab2016b

% Copyright © 2018-2020 Tomas Lukes, Kristin Grussmayer, Adrien Descloux
% École Polytechnique Fédérale de Lausanne,
% Laboratory of Nanoscale Biology, http://lben.epfl.ch/
% kristin.grussmayer@epfl.ch, adrien.descloux@epfl.ch, tomas.lukes@epfl.ch
 
% This file is part of multicolorSOFI.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%
%%% Initialization
%

addpath(genpath('funcs'));
config;

%% Load data
disp('Loading reflected data')
imstack1 = load_tifFile([settings.io.pn, filesep, settings.io.fn1, ...
    settings.io.fileExtension],settings.sys.sub);
disp('Loading transmitted data')
imstack2 = load_tifFile([settings.io.pn, filesep, settings.io.fn2, ...
    settings.io.fileExtension],settings.sys.sub);
% flip experimental data from second camera
if strcmp(settings.io.dataType, 'exper')
    imstack2=flip(imstack2,2);
end

%
%%% Main Program
%

%% Run multicolor SOFI
[umx_c2, im_dec, im_lin, mean1c, mean2tc] = multicolor_sofi_interface(...
    imstack1, imstack2, settings);

%
%%% Plotting, visualisation
%

%% Show linearized deconvolved and unmixed images as RGB images
im_rgb_sofi2 = mergeToRgb(mean(umx_c2(:,:,:,3),3), ...
    mean(umx_c2(:,:,:,2),3), mean(umx_c2(:,:,:,1),3));

im_rgb_sofi2_lin = mergeToRgb(mean(im_lin(:,:,:,3),3), ...
    mean(im_lin(:,:,:,2),3), mean(im_lin(:,:,:,1),3));

figure(5)
subplot(131);
imshowpair(mean1c,mean2tc,'falsecolor');axis equal;axis tight;
title('Mean overlay')
subplot(132);
imshow(im_rgb_sofi2,[])
title('SOFI 2 unmixed')
subplot(133);
imshow(im_rgb_sofi2_lin,[]);
title('SOFI 2 lin dec unmixed') %unmixed, deconvolved and linearised
set(gcf,'position',[465   555   1200   350]);
if settings.io.figshow == 0; set(gcf,'visible', 'off'); end

r_lin = im_rgb_sofi2_lin(:,:,1);
g_lin = im_rgb_sofi2_lin(:,:,2);
b_lin = im_rgb_sofi2_lin(:,:,3);

figure(6)
subplot(131);
rf = zeros(size(r_lin,1),size(r_lin,2),3); 
rf(:,:,1) = r_lin;
imshow(rf./max(rf(:)));
title('Red channel');
subplot(132);
gf = zeros(size(g_lin,1),size(g_lin,2),3); 
gf(:,:,2) = g_lin;
imshow(gf./max(gf(:)));
title('Green channel');
subplot(133);
bf = zeros(size(b_lin,1),size(b_lin,2),3); 
bf(:,:,3) = b_lin;
imshow(bf./max(bf(:)));
title('Blue channel');
set(gcf,'position',[465   555   1200   350])
if settings.io.figshow == 0; set(gcf,'visible', 'off'); end

%% TEST expected results

expected = load('test_files/simulated_data_umx_c2.mat');
if isequal(expected.umx_c2_all, umx_c2) == 1
    disp('Unmixed cumulants TEST PASSED')
else
    disp('Unmixed cumulants TEST FAILED')
end
expected = load('test_files/simulated_data_im_dec.mat');
if isequal(expected.im_dec, im_dec) == 1
    disp('Image deconvolution TEST PASSED')
else
    disp('Image deconvolution TEST FAILED')
end
expected = load('test_files/simulated_data_im_lin.mat');
if isequal(expected.im_lin, im_lin) == 1
    disp('Image linearisation TEST PASSED')
else
    disp('Image linearisation TEST FAILED')
end
expected = load('test_files/simulated_data_rgb_sofi.mat');
if isequal(expected.im_rgb_sofi2_lin, im_rgb_sofi2_lin) == 1
    disp('Final RGB image equality TEST PASSED')
else
    disp('Final RGB image equality TEST FAILED')
end

%% Saving
settings.io.outputpath = pwd;
if settings.io.figsave == 1 
    disp('Saving results')
    results_folder = ['results_',datestr(now,'yyyy-mm-dd_HHMMSS')];
    mkdir(settings.io.outputpath,results_folder);
    outputfn = regexprep([settings.io.fn1],'_reflected','','ignorecase');

    writeTIFF(mean1c./max(mean1c(:)),[settings.io.outputpath,...
        filesep, results_folder,filesep,outputfn,...
        '_mean_reflected_channel'])
    writeTIFF(mean2tc./max(mean2tc(:)),[settings.io.outputpath,...
        filesep,results_folder,filesep,outputfn,...
        '_mean_transmitted_channel'])
        
    writeRGBTIFF(uint8(255.*im_rgb_sofi2),[settings.io.outputpath,...
        filesep,results_folder,filesep,outputfn,...
        '_sofi2_unmixed_rgb_composite'])
    writeRGBTIFF(uint8(255.*im_rgb_sofi2_lin),[settings.io.outputpath,...
        filesep,results_folder,filesep,outputfn,...
        '_sofi2_un  mixed_linearised_rgb_composite'])
    
    writeTIFF(squeeze(mean(umx_c2,3)),[settings.io.outputpath,...
        filesep,results_folder,filesep,outputfn,'_sofi2_unmixed'])
    writeTIFF(squeeze(mean(im_lin,3)),[settings.io.outputpath,...
        filesep,results_folder,filesep,outputfn,'_sofi2_unmixed_linearised'])

    saveFigure(5,[settings.io.outputpath,filesep,results_folder],...
        [outputfn,'_results_rgb_composite'],'tif')
    saveFigure(6,[settings.io.outputpath,filesep,results_folder],...
        [outputfn,'_results_channels'],'tif')
    saveSettingsTxt(settings,[settings.io.outputpath,filesep,...
        results_folder,filesep,outputfn,'_settings.txt'])
end
