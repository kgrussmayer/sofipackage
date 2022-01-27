function results = estimateton(stack,settings,results)
%  estimate ton by calculating 2nd order cumulant images with different
%  time shifts and then meantrace over these images
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings
%           results ... struct for passing intermediate results across the
%           processing chain
% output:   results ... struct ton estimate and fitted meantrace

tau=unique(round(logspace(0,log10(settings.ton.subseqlength/2),settings.ton.numtau)));

for nt=1:numel(tau)

    stack3Dt = cat(4,stack(:,:,(1: settings.ton.subseqlength)),...
                 stack(:,:, 1+tau(nt): settings.ton.subseqlength+tau(nt))); 

    stack3Dt = permute(stack3Dt,[1 2 4 3]);
    [sofi,grid]=sofiCumulants2Dt(stack3Dt,1,[],[],1:2);
    soficor(nt)=mean2(abs(sofi{2}(100:150,100:150,2)));
    sofitau(:,:,nt) = sofi{2}(:,:,2);
end

x=(0:length(soficor)-1);
mitrace=soficor;

if contains(struct2array(ver), 'Curve Fitting Toolbox')
    % this version uses the curve fitting toolbox
        s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[Inf,Inf,Inf],...
               'Startpoint',[1 1 1]);
        f = fittype('a*exp(-x/b)+c','options',s);

    [c2,gof2] = fit(x',mitrace(1:end)',f);
    a = c2.a;
    b = c2.b;
    c = c2.c;
    
    fitcurve = a*exp(-x/b)+c;

    results.ton.mitrace = mitrace;
    results.ton.fitinfo = c2;
    results.ton.fitcurve = fitcurve;
    results.ton.xaxis = x;
elseif contains(struct2array(ver), 'Statistics and Machine Learning Toolbox')
    % this version uses fitnlm from Statistics package 
    % Convert X and Y into a table
    tbl = table(x', mitrace');
    % define the model as Y = a*exp(-t/b)+c
    modelfun = @(b,t)  b(1)*exp(-t(:, 1)/b(2))+ b(3) ;  
    beta0 = [1, 1, 1];  
    % actual model computation
    mdl = fitnlm(tbl, modelfun, beta0);

    % Extract the coefficient values from the the model object.
    % The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
    coefficients = mdl.Coefficients{:, 'Estimate'};

    a = coefficients(1);
    b = coefficients(2);
    c = coefficients(3);

    fitcurve = a*exp(-x/b)+c;
    
    results.ton.mitrace = mitrace;
    results.ton.fitinfo = mdl;
    results.ton.fitcurve = fitcurve;
    results.ton.xaxis = x;

    disp('utils.bleachcorr.estimateton: Fit using fitnlm without parameter bounds: Please check the fit quality.')
else
    disp("Fitting without packages to be added")
end
    
% plot and save meantrace figure
if settings.io.figs ==1   
    figure;
    plot(x',mitrace);hold on;plot(x,fitcurve,'k');
    ylabel('Decorelation curve');xlabel('Number of frames');
    title(['Ton estimation, average Ton ',num2str(b),' frames']); 
    
    if settings.io.figsave ==1
    saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
        [settings.io.imageName,'_ton'],settings.io.figformat,0);
    end
end
% eof