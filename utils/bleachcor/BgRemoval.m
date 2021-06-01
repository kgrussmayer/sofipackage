function [DataCorrected,Background] = BgRemoval(Data)

% *************************************************************************
% Interface to remove background of a single-molecule fluorescence 
% microscopy imaging data set. The method is based on a filtering of the 
% Principal Component Analysis (PCA) loadings. SOFI can then be applied 
% onto the corrected data.
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
Continue = 1;
comp = 1;
Stop = 'n';
[m,n,o] = size(Data);
DataMatrix = reshape(Data,m*n,o);
screensize = get(groot,'Screensize');

% Error messages when providing wrong input.
if o == 1
    error('Please provide the input data in its data cube form (x, y, time)');
end

% Perform calculations until user specifies to stop.
while Continue;
    
    % Visualization of the calculations.
    reverseString = '';
    msg = sprintf(['Calculating background contribution ',num2str(comp),'. Please wait a moment...']);
    fprintf([reverseString, msg]);
    reverseString = repmat(sprintf('\b'), 1, length(msg));
    
    % Calculate the background component of the data.
    [~,~,V] = svds(DataMatrix,1);
    
    % Reset the automatic selection.
    Smoothing = 1;
    Autom = 0;
    
    % Continue with the same component until user specifies to stop.
    while Smoothing == 1
        
        % Matrix used for visualization of the results.
        DataMatrix2 = DataMatrix-min(DataMatrix,1);
        
        % Plot the results.
        % The first figure is the PCA loadings (blue line) to which the
        %   filtered component will be added later (green line).
        % The second figure is the reconstructed data.
        % The third figure is the mean trace of the reconstructed data.
        %   This can be used to decide when to stop (fairly constant line
        %   with blinking means that the background has been removed).
        h = figure(1);
        set(h,'Position',[screensize(3)*0.25 screensize(4)*0.15 screensize(3)*0.5 screensize(4)*0.7],'NumberTitle','off','Name','Background removal by PCA: analysis');
        subplot(4,4,[1:8]);h1=plot(V,'b');hold on;title(['Background contribution ',num2str(comp)]);
        subplot(4,4,[9 10 13 14]);imagesc(mean(reshape(DataMatrix2,m,n,o),3));colorbar;
        if comp == 1
            title('Raw Data');
        else
            title('Mean image of corrected data');
        end
        if m == n
            axis square;
        end
        subplot(4,4,[11 12 15 16]);plot(mean(DataMatrix2,1));axis square
        if comp == 1
            title('Mean time trace of raw data');
        else
            title('Mean time trace of corrected data');
        end
        shg
        
        % Perform the automatic parameter selection in case it has not been
        % done before. This is done for every background component.
        if Autom == 0
            
            % Ask to stop the analysis, based on looking at the background
            % component.
            if comp ~= 1
                
                fprintf(reverseString);
                Stop = input('Stop analysis (y/n)? ','s');
                
                % Visualization of the calculations.
                msg = sprintf(['\nStop analysis (y/n)? ']);
                reverseString = repmat(sprintf('\b'), 1, length(msg)+1);
                
            end
            
            % Exit the loop and analysis will be stopped.
            if Stop == 'y'
                break
            end
            
            % Visualization of the calculations.
            msg = sprintf(['Proposing an optimal smoothing parameter (Lambda). Please wait a moment...']);
            fprintf([reverseString, msg]);
            reverseString = repmat(sprintf('\b'), 1, length(msg));
                        
            % Perform the automatic parameter selection (100 values are
            % checked) by using a leave-one-out cross-validation method.
            % Only part of the signal is taken, as the blinking can be
            % correlated, which would induce errors in the parameter
            % selection.
            Lambda = logspace(-8,10,100);
            for i = 1:size(Lambda,2)
                [~,CV(i)] = SmoothPCA(V(1:5:end), Lambda(i));
            end
            [~, ind] = min(CV);
            
            % Correction for thinning of the signal.
            %   (Lambda_new = Lambda_old * f^(2*d)
            % f is the thinning factor (5) and d is fixed to be 2.
            Lambda = Lambda(ind)*(5)^(2*2);
            
            % Limitation of Lambda, to avoid instabilities in the solution.
            if Lambda > 1e12
                Lambda = 1e12;
            end
            
            % Calculation of the filtered (smoothed) background component.
            VSm(:,comp) = SmoothPCA(V, Lambda);
        
        % If automatic selection has been done before, but results were not
        % satisfactory, manuel filtering can be done.
        else
            
            % Manual filtering Lambda input
            fprintf(reverseString);
            LambdaInput = input(['Input 0 to stop the analysis or give a smoothing parameter (previous value: ',num2str(round(Lambda)),'): '],'s');
            
            % Visualization of the calculations.
            msg = sprintf(['\nInput 0 to stop the analysis or give a smoothing parameter (previous value: ',num2str(round(Lambda)),'): ']);
            reverseString = repmat(sprintf('\b'), 1, length(msg)+length(LambdaInput));
            
            % Convert the input to an integer.
            Lambda = str2num(LambdaInput);
            
            % Limitation of Lambda, to avoid instabilities in the solution.
            if Lambda > 1e12
                msg = sprintf(['Attention: smoothing parameter (Lambda) has been limited to 10^12 to avoid instabilities.']);
                fprintf([reverseString, msg]);
                reverseString = repmat(sprintf('\b'), 1, length(msg));
                
                Lambda = 1e12;
            end
            
            % Calculate the manual filtering
            VSm(:,comp) = SmoothPCA(V, Lambda);
            
        end
        
        % If analysis should be stopped or manual filtering has a 0 value
        % for the smoothing parameter, the loop will be exited and analysis
        % will stop.
        if Stop == 'y' || Lambda == 0
            break
        end
        
        % Recalculation of the matrix to show the visual result of the
        % filtering of the background component.
        DataMatrix2 = (DataMatrix - DataMatrix*VSm(:,comp)*VSm(:,comp)')-min(DataMatrix - DataMatrix*VSm(:,comp)*VSm(:,comp)',1);
        
        % Update the figure.
        %   Add the filtered background component to the first figure.
        %   Update the reconstructed data in the second figure.
        %   Update the mean trace of the reconstructed data in the third
        %   figure.
        figure(1);
        subplot(4,4,[1:8]);h2=plot(VSm(:,comp),'g','LineWidth',2);hold off
        subplot(4,4,[9 10 13 14]);imagesc(mean(reshape(DataMatrix2,m,n,o),3));colorbar;title('Mean image of corrected data');
        if m == n
            axis square;
        end
        subplot(4,4,[11 12 15 16]);plot(mean(DataMatrix2,1));axis tight;axis square;title('Mean time trace of corrected data');
        shg
        
        % If automatic parameter selection has been performed, ask if
        % results are satisfactory or not.
        if Autom == 0
            
            fprintf(reverseString);
            Smoothing = input(['Accept automatic smoothing (Lambda = ',num2str(round(Lambda)),') and go to next background contribution (0) / Continue filtering current background contribution (1): ']);
            Autom = 1;
            
            % Visualization of the calculations.
            msg = sprintf(['\nAccept automatic smoothing (Lambda = ',num2str(round(Lambda)),') and go to next background contribution (0) / Continue filtering current background contribution (1): ']);
            reverseString = repmat(sprintf('\b'), 1, length(msg)+1);
        
        % If no automatic parameter selection was performed, ask if manual 
        % selection is satisfactory.
        else
            
            fprintf(reverseString);
            Smoothing = input('Go to next background contribution (0) / Continue filtering current background contribution (1): ');
            
            % Visualization of the calculations.
            msg = sprintf(['\nGo to next background contribution (0) / Continue filtering current background contribution (1): ']);
            reverseString = repmat(sprintf('\b'), 1, length(msg)+1);
        
        end
        
    end
    
    % The loop will be exited and analysis will stop.
    if Stop == 'y' || Lambda == 0
       break
    end
    
    % Visualization of the calculations.
    msg = sprintf(['Finished treating background contribution ',num2str(comp),' with a smoothing parameter (Lambda) of ',num2str(round(Lambda)),'.']);
    fprintf([reverseString, msg,'\n']);
    
    % Recalculate the reconstructed data, after the filtering of each
    % background component.
    DataMatrix = DataMatrix - DataMatrix*VSm(:,comp)*VSm(:,comp)';
    
    comp = comp + 1;
    
end

% Visualization of the calculations.
if comp-1 == 1
    msg = sprintf(['Finished background removal procedure. A total of ',num2str(comp-1),' background contribution was removed.']);
else
    msg = sprintf(['Finished background removal procedure. A total of ',num2str(comp-1),' background contributions were removed.']);
end
fprintf([reverseString, msg,'\n']);

% When analysis is finished, reconstruct the final data and background.
DataCorrected = reshape(DataMatrix-min(DataMatrix,1),m,n,o);
Background = reshape((reshape(Data,m*n,o)-DataMatrix)-min((reshape(Data,m*n,o)-DataMatrix),1),m,n,o);

% Show the final results of the background removal.
%   The first figure shows the raw data.
%   The second figure shows the corrected data.
set(h,'Position',[screensize(3)*0.25 screensize(4)*0.15 screensize(3)*0.5 screensize(4)*0.7],'NumberTitle','off','Name','Background removal by PCA: final results');
subplot(1,2,1);imagesc(mean(Data,3));title('Raw data')
if m == n
	axis square;
end
subplot(1,2,2);imagesc(mean(DataCorrected,3));title('Corrected data')
if m == n
	axis square;
end
if m ~= n
    truesize(h);
end
shg
    