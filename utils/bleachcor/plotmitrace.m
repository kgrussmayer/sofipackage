function [settings,results] = plotmitrace(stack,settings)
% calculate meantrace and estimate average bleaching lifetime 
% exponential fitting to the curve given by the decay of average intensity for each frame 
% depending on the toolbox availability, different fitting options are used
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings

frames = size(stack,3);
mitrace=squeeze(mean(mean(stack,1),2));
mitrace = mitrace./max(mitrace);

if contains(struct2array(ver), 'Curve Fitting Toolbox')
    % this version uses the curve fitting toolbox
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,Inf,Inf],...
                   'Startpoint',[1 1 1]);
    f = fittype('a*exp(-x/b)+c','options',s);

    [c2,~] = fit((0:frames-1)',mitrace(1:end),f);
    a = c2.a;
    b = c2.b;
    c = c2.c;

    x=(0:frames-1);
    fitcurve = a*exp(-x/b)+c;  
    results.blcor.fitinfo = c2;
elseif contains(struct2array(ver), 'Statistics and Machine Learning Toolbox')
    % this version uses fitnlm from Statistics package 
    x=(0:frames-1);

    % Convert X and Y into a table
    tbl = table(x', mitrace);
    % define the model as Y = a*exp(-t/b)+c
    modelfun = @(b,t)  b(1)*exp(-t(:, 1)/b(2))+ b(3) ;  
    beta0 = [1, 100, 1];  
    % actual model computation
    mdl = fitnlm(tbl, modelfun, beta0);

    % Extract the coefficient values from the the model object.
    % The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
    coefficients = mdl.Coefficients{:, 'Estimate'};

    a = coefficients(1);
    b = coefficients(2);
    c = coefficients(3);

    fitcurve = a*exp(-x/b)+c;
    results.blcor.fitinfo = mdl;
    disp('utils.bleachcorr.plotmitrace: Fit using fitnlm without parameter bounds: Please check the fit quality.')
else
    % without any toolbox
    % fit inputs in the same shape: frame index as time, mean intensity
    x = (0:frames-1);
    y = mitrace';
    % exponential function with offset a*exp(-l*t)+c
    % with parameters vector p = [a, l, c]
    f = @(p, t) p(1).*exp(-p(2).*t)+p(3);
    % objective function to minimize: ordinary least squares for now
    cost = @(p) sum((f(p, x) - y).^2);
    opts = optimset('MaxFunEvals', 50000, 'MaxIter', 10000);
    p0 = [1 1/100 1];
    P = fminsearch(cost, p0, opts);
    % assign estimated parameters
    a = P(1);
    b = 1/P(2);
    c = P(3);
    % evaluate fitted function
    fitcurve = f(P, x);
    results.blcor.fitinfo = [];
    disp('utils.bleachcorr.plotmitrace: Fit using fminsearch without parameter bounds: Please check the fit quality.')
end

if settings.blcor.MaxCorrSamp > size(stack,3)
    settings.blcor.MaxCorrSamp = size(stack,3);
end
[~,corelf] = stackcorel(stack,settings.blcor.MaxCorrSamp);

results.blcor.mitrace = mitrace;
results.blcor.fitcurve = fitcurve;
results.blcor.xaxis = x;
results.blcor.corelf = corelf;

 
% plot and save meantrace figure
if settings.io.figs ==1 
    figure;plot(x',mitrace);hold on;plot(x,fitcurve,'k');
    ylabel('Mean intensity');xlabel('Number of frames');
    title(['Bleaching estimation, average bleaching lifetime ',num2str(b),' frames']); 

    if settings.io.figsave ==1 
    saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
        [settings.io.imageName,'_meantrace'],settings.io.figformat,0);
    end
end
% eof