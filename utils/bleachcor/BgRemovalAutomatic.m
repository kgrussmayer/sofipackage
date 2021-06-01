function [DataCorrected,Background] = BgRemovalAutomatic(Data,ncomp)

% *************************************************************************
% Fully automatic routing to remove background of a single-molecule 
% fluorescence microscopy imaging data set. The method is based on a 
% filtering of the Principal Component Analysis (PCA) loadings. SOFI can 
% then be applied onto the corrected data.
%
% *************************************************************************
% Authors:  
%       - HUGELIER, S.      (1)
%       - VITALE, R.        (1,2)
%       - RUCKEBUSCH, C.    (1)
%           
% (1):  LAboratoire de Spectrochimie Infrarouge et Raman (LASIR)
%       Université de Lille 1, UMR CNRS 8516
%       Bât C5, Cité Scientifique
%       59655 Villeneuve d'Ascq - France
% (2):  Laboratory of Photochemistry and Spectroscopy
%       Katholieke Universiteit Leuven
%       B-3001 Heverlee, Belgium
%
% *************************************************************************
% Input:
%   Data:   Three-dimensional fluorescence image (x,y,time).
%   ncomp:  The number of expected background contributions.
%
% Output:
%   DataCorrected:  Corrected data (to which SOFI can be applied).
%   Background:   	Removed background.
%
% *************************************************************************

% Close the figure needed during the analysis.
if any(findall(0,'Type','Figure')==1)
    close(1);
end

% clear the command window
clc

% Initialization for the calculations.
comp = 1;
[m,n,o] = size(Data);
DataMatrix = reshape(Data,m*n,o);
screensize = get(groot,'Screensize');
if nargin < 2
    ncomp = 4;
end

% Error messages when providing wrong input.
if o == 1
    error('Please provide the input data in its data cube form (x, y, time)');
end
if ncomp > 8
    error('Only a maximum of 8 expected background contributions is allowed.');
end

% Perform calculations for the specified number of background 
% contributions.
while comp <= ncomp;
    
    % Visualization of the calculations.
    reverseString = '';
    msg = sprintf(['Calculating background contribution ',num2str(comp),'. Please wait a moment...']);
    fprintf([reverseString, msg]);
    reverseString = repmat(sprintf('\b'), 1, length(msg));
    
    % Calculate the background component of the data.
    [~,~,V(:,comp)] = svds(DataMatrix,1);
    
    % Visualization of the calculations.
    msg = sprintf(['Calculating the optimal smoothing parameter (Lambda) for background contribution ',num2str(comp),'. Please wait a moment...']);
    fprintf([reverseString, msg]);
    reverseString = repmat(sprintf('\b'), 1, length(msg));
        
    % Perform the automatic parameter selection (100 values are checked) by 
    % using a leave-one-out cross-validation method. Only part of the 
    % signal is taken, as the blinking can be correlated, which would 
    % induce errors in the parameter selection.
    Lambda = logspace(-8,10,100);
    for i = 1:size(Lambda,2)
    	[~,CV(i)] = SmoothPCA(V(1:5:end,comp), Lambda(i));
    end
    [~, ind] = min(CV);
            
    % Correction for thinning of the signal.
    %   (Lambda_new = Lambda_old * f^(2*d)
    %   f is the thinning factor (5) and d is fixed to be 2.
    Lambda = Lambda(ind)*(5)^(2*2);
            
    % Limitation of Lambda, to avoid instabilities in the solution.
    if Lambda > 1e12
    	Lambda = 1e12;
    end
            
    % Calculation of the filtered (smoothed) background component.
    VSm(:,comp) = SmoothPCA(V(:,comp), Lambda);
    
    % Visualization of the calculations.
    msg = sprintf(['Finished treating background contribution ',num2str(comp),' with a smoothing parameter (Lambda) of ',num2str(round(Lambda)),'.']);
    fprintf([reverseString, msg,'\n']);
    
    % Recalculate the reconstructed data, after the filtering of each
    % background component.
    DataMatrix = DataMatrix - DataMatrix*VSm(:,comp)*VSm(:,comp)';
    
    comp = comp + 1;
    
end

% Visualization of the calculations.
msg = sprintf(['Finished background removal procedure. A total of ',num2str(comp-1),' background contributions were removed.']);
fprintf([msg,'\n']);

% When analysis is finished, reconstruct the final data and background.
DataCorrected = reshape(DataMatrix-min(DataMatrix,1),m,n,o);
Background = reshape((reshape(Data,m*n,o)-DataMatrix)-min((reshape(Data,m*n,o)-DataMatrix),1),m,n,o);

% Show the final results of the background removal.
%   The first set of figures shows the background contributions and the
%   filtered parts.
%   The second figure shows the raw data.
%   The third figure shows the corrected data.
h = figure(1);
set(h,'Position',[screensize(3)*0.25 screensize(4)*0.15 screensize(3)*0.5 screensize(4)*0.7],'NumberTitle','off','Name','Background removal by PCA: final results');
for i = 1:ncomp
    if ncomp <= 4
        subplot(4,4,[(i-1)*2+1:i*2]);plot(V(:,i),'b');hold on;plot(VSm(:,i),'g','LineWidth',2);title(['Background contribution ',num2str(i)]);
    else
        subplot(4,4,i);plot(V(:,i),'b');hold on;plot(VSm(:,i),'g','LineWidth',2);title(['Background contribution ',num2str(i)]);
    end
end
subplot(4,4,[9 10 13 14]);imagesc(mean(Data,3));title('Raw data')
if m == n
	axis square;
end
subplot(4,4,[11 12 15 16]);imagesc(mean(DataCorrected,3));title('Corrected data')
if m == n
	axis square;
end
if m ~= n
    truesize(h);
end
shg
    