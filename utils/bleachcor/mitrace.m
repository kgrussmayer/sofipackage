function settings = mitrace(stack,settings)
% calculate meantrace and estimate average bleaching lifetime 
% gaussian fitting to the curve given by the decay of average intensity for each frame 
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings

frames = size(stack,3);

mitrace=squeeze(mean(mean(stack,1),2));
mitrace = mitrace./max(mitrace);
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

settings.blcor.fitinfo = c2;
settings.blcor.fitcurve = fitcurve;

% plot and save meantrace figure
figure,
plot(x',mitrace);hold on;plot(x,fitcurve,'k');
ylabel('Mean intensity');xlabel('Number of frames');
title(['Bleaching estimation, average bleaching lifetime ',num2str(b),' frames']); 

saveFigure(gcf,[settings.io.outputpath,filesep,'figs'],[settings.io.imageName,'_meantrace'],settings.io.figformat,0);
end
% eof