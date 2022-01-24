function [settings,results] = plotmitrace(stack,settings)
% calculate meantrace and estimate average bleaching lifetime 
% gaussian fitting to the curve given by the decay of average intensity for each frame 
% this version uses fitnlm from Statistics package and was written by
% Monika Paw³owska, based on https://se.mathworks.com/matlabcentral/answers/91159-how-do-i-fit-an-exponential-curve-to-my-data
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings

frames = size(stack,3);
mitrace=squeeze(mean(mean(stack,1),2));
mitrace = mitrace./max(mitrace);

X=(0:frames-1);

% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X', mitrace);
% Define the model as Y = a*exp(-x/b)+c
modelfun = @(b,x)  b(1)*exp(-x(:, 1)/b(2))+ b(3) ;  
beta0 = [1, 100, 1];  
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);

% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'};

a = coefficients(1);
b = coefficients(2);
c = coefficients(3);

fitcurve = a*exp(-X/b)+c;

if settings.blcor.MaxCorrSamp > size(stack,3)
    settings.blcor.MaxCorrSamp = size(stack,3);
end
[~,corelf] = stackcorel(stack,settings.blcor.MaxCorrSamp);

results.blcor.mitrace = mitrace;
results.blcor.fitinfo = mdl;
results.blcor.fitcurve = fitcurve;
results.blcor.xaxis = X;
results.blcor.corelf = corelf;
 
% plot and save meantrace figure
if settings.io.figs ==1 
    figure;plot(X',mitrace);hold on;plot(X,fitcurve,'k');
    ylabel('Mean intensity');xlabel('Number of frames');
    title(['Bleaching estimation, average bleaching lifetime ',num2str(b),' frames']); 

    if settings.io.figsave ==1 
    saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
        [settings.io.imageName,'_meantrace'],settings.io.figformat,0);
    end
end
% eof