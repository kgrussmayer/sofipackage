%%% batch code for MCSOFI processing to determine the splitting across the 
%%% dichroic from experimental 1 color labeled data; calculates the reflection and
%%% transmission coefficient
%%% kristin.grussmayer@epfl.ch 10/2018

close all; clear; 

%%% Initialization

addpath(genpath('..\utils'));
addpath(genpath('funcs'));
config_mc;

% select folder to batch process
% this folder should contain folder with beads calibration data for the
% dichroic and pairs of images for the reflection and transmission channel

%ATTENTION: CURRENTLY rlist AND tlist EXCHANGED IN MAIN LOOP DUE TO WRONG
%FILE NAMES!
% pathname to extra folder for pairs of files
pname = 'G:\Kristin\20181017\Alexa568';
settings.io.pn = pname;
% pathname to calibration folder
isCal = 1;
settings.io.pnc = 'G:\Kristin\20181017\calibration-ZT594';

list = dir(pname);
list(1:2) = []; % remove . and .. path from the list

% delete files that don't need processing
for k = length(list):-1:1 % loop backwards because we are changing the list on the fly
    if strfind(list(k).name,'results')
        list(k) = []; % remove results folder from the list
    elseif isempty(strfind(list(k).name,'.tiff'))&& isempty(strfind(list(k).name,'.tif'))
        list(k) = []; % remove extra folders from tomas procesing
    elseif strfind(list(k).name,'.mat')
        list(k) = [];
    end
end
if isCal == 0
    disp('Calibration file not found !');
    disp('Please make sure the folder is properly organized!'); 
else    

% make list of paired files
ex = length(list(1).name)-strfind(list(1).name,'.');
idMax = 0; listC = [];
for k = 1:length(list)
    temp = str2double(list(k).name(end-ex-3:end-ex-1));
    if temp > idMax
        idMax = temp;
    end
end

rList = []; tList = [];
for k = 0:idMax 
    id = num2str(k); while length(id) < 3; id = ['0',id];end
    id = ['_',id];
    for h = length(list):-1:1
       if strfind(list(h).name,id) % if we find proper ID number
           if strfind(list(h).name,'reflected')
               rList(end+1).name = list(h).name;
               rList(end).id = id(2:end);
           elseif strfind(list(h).name,'transmitted')
               tList(end+1).name = list(h).name;
               tList(end).id = id(2:end);
           end
       end
    end
end

if length(tList) ~= length(rList)
    warning('Error number of reflected and transmitted files are different');
end

% summary of folder to be processed
disp(['Selected folder : ',settings.io.pn])
disp(['Calibration folder : ',settings.io.pnc])
disp(['Number of files detected : ',num2str(length(rList))])
end

%% Load the calibration files
imc1 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc1, settings.io.fileExtension]);
imc2 = load_tifFile([settings.io.pnc, filesep, settings.io.fnc2, settings.io.fileExtension]);
imc2 = flip(imc2,2);

% imc1(imc1==0)=min(imc1(imc1>0));
% imc2(imc2==0)=min(imc2(imc2>0));

% what is indcom used for?
[~,indcom]=max(squeeze(mean(max(imc1))).*squeeze(mean(max(imc2))));
im_fix=imc1(:,:,indcom);
im_mov=imc2(:,:,indcom);

show2subs(im_fix,im_mov,1)
close(1)

%%% Run calibration (register beads images from 2 color channels)
disp('Compute calibration')
settings.cal.roix = []; settings.cal.roiy = [];
cal = runcalib_mcsofi(imc1,imc2,settings.sys,settings.cal,settings.io); 

%% MAIN PROCESSING LOOP

R = zeros(length(rList), 1);

for k = 1:length(rList) 
    disp(['Processing file # ',num2str(k),', ID : ',rList(k).id])

%% Load stack of images to be transformed (1st stack) 

settings.io.fn1 = tList(k).name;
settings.io.fn2 = rList(k).name;
disp('Loading reflected data')
imstack1 = load_tifFile([settings.io.pn, filesep, settings.io.fn1],settings.sys.sub);
disp('Loading transmitted data')
imstack2 = load_tifFile([settings.io.pn, filesep, settings.io.fn2],settings.sys.sub);
imstack2=flip(imstack2,2);

%% Transform stack of images
Nframes = size(imstack2,3);
disp('Transform of imstack2')
% apply transform to data
for ii = 1:Nframes
    im_mov = imstack2(:,:,ii);    
    moving_reg = imwarp(im_mov, cal.tf{1},'OutputView',imref2d(size(imstack1)));
    imstack2t(:,:,ii) = moving_reg;
end

settings.cal.shiftx = -cal.tf{1}.A(1);
settings.cal.shifty = -cal.tf{1}.B(1);

% compute coregistration mask and proper crop
mask = not(imstack2t(:,:,1) == 0);
[r,c] = find(mask>0);
rect = [min(c)+1 min(r)+1 max(c)-min(c)-2 max(r)-min(r)-2];

% SOFI processing requires even dimensions or you gets weird cropping
if mod(rect(3),2) == 0; rect(3) = rect(3)-1; end
if mod(rect(4),2) == 0; rect(4) = rect(4)-1; end

settings.cal.roix = [rect(1) rect(3)];
settings.cal.roiy = [rect(2) rect(4)];

% crop image stack 
imstack1crop = imstack1(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:);
% crop image stack 
imstack2crop = single(imstack2t(rect(2):rect(2) + rect(4),rect(1):rect(1)+rect(3),:));

% display rgb image to check coregistration on data
fake = [];
fake(:,:,2) = 255.*imstack1crop./max(imstack1crop(:)); fake(:,:,3) = 255.*imstack2crop./max(imstack2crop(:)); fake(:,:,1) = 255.*imstack2crop./max(imstack2crop(:));
figure(2);imagesc(uint8(fake))

%% calculate the R/T coefficient
% mean pixel background in reflected & transmitted camera frames with laser shutters shut
% measurements mean of 10 different movies
rMeanPixelBkg = 102.25;
tMeanPixelBkg = 100.65;
% calculate the percentage of reflected light
% signal = sum(pixel intensities) - #pixels * meanPixelBkg
[px,py] = size(imstack1crop);
test = im2double(imstack1crop);
signal_reflected = sum(sum(imstack1crop))- px*py*rMeanPixelBkg;
signal_transmitted = sum(sum(imstack2crop))- px*py*tMeanPixelBkg;
r = signal_reflected/ (signal_reflected + signal_transmitted);
R(k, 1) = r;
end
dlmwrite(fullfile(pname, 'R_mean_std'), [R; mean(R, 1); std(R, 1)]');
T = 1-R;
dlmwrite(fullfile(pname, 'T_mean_std'), [T; mean(T, 1); std(T, 1)]');


