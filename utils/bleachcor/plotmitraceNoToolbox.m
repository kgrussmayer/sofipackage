function [settings,results] = plotmitrace(stack,settings)
% calculate meantrace and estimate average bleaching lifetime 
% gaussian fitting to the curve given by the decay of average intensity for each frame 
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings

frames = size(stack,3);
mitrace = squeeze(mean(mean(stack,1),2));
mitrace = mitrace./max(mitrace);

% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,0],...
%                'Upper',[Inf,Inf,Inf],...
%                'Startpoint',[1 1 1]);
% f = fittype('a*exp(-x/b)+c','options',s);
% 
% [c2,~] = fit((0:frames-1)',mitrace(1:end),f);
% a = c2.a;
% b = c2.b;
% c = c2.c;
%
% x = (0:frames-1);
% fitcurve = a*exp(-x/b)+c

% fit inputs in the same shape: frame index as time, mean intensity
x = (0:frames-1);
y = mitrace';
% exponential function with offset a*exp(-l*t)+c
% with parameters vector p = [a, l, c]
f = @(p, t) p(1).*exp(-p(2).*t)+p(3);
% objective function to minimize: ordinary least squares for now
cost = @(p) sum((f(p, x) - y).^2);
opts = optimset('MaxFunEvals', 50000, 'MaxIter', 10000);
p0 = [1 1 1];
P = fminsearch(cost, p0, opts)
% evaluate fitted function
fitcurve = f(P, x);

if settings.blcor.MaxCorrSamp > size(stack,3)
    settings.blcor.MaxCorrSamp = size(stack,3);
end
[~,corelf] = stackcorel(stack,settings.blcor.MaxCorrSamp);

results.blcor.mitrace = mitrace;
results.blcor.fitinfo = [];
results.blcor.fitcurve = fitcurve;
results.blcor.xaxis = x;
results.blcor.corelf = corelf;
 
% plot and save meantrace figure
if settings.io.figs ==1 
    figure;plot(x',mitrace);hold on;plot(x,fitcurve,'k');
    ylabel('Mean intensity');xlabel('Number of frames');
    title(['Bleaching estimation, average bleaching lifetime ', num2str(1/P(2)),' frames']); 

    if settings.io.figsave ==1 
    saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],...
        [settings.io.imageName,'_meantrace'],settings.io.figformat,0);
    end
end
% eof