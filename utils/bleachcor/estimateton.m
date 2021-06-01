function results = estimateton(stack,settings,results)
%  estimate ton by calculating 2nd order cumulant images with different
%  time shifts and then meantrace over these images
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings
%           results ... struct for passing intermediate results across the
%           processing chain
% output:   

tau=unique(round(logspace(0,log10(settings.ton.wsize/2),settings.ton.numtau)));

for nt=1:numel(tau)

stack3Dt = cat(4,stack(:,:,(1: settings.ton.wsize)),...
                 stack(:,:, 1+tau(nt): settings.ton.wsize+tau(nt))); 

stack3Dt = permute(stack3Dt,[1 2 4 3]);
[sofi,grid]=sofiCumulants2Dt(stack3Dt,1,[],[],1:2);
soficor(nt)=mean2(abs(sofi{2}(100:150,100:150,2)));
sofitau(:,:,nt) = sofi{2}(:,:,2);
end

mitrace=soficor;
%         mitrace = mitrace./max(mitrace);
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[Inf,Inf,Inf],...
               'Startpoint',[1 1 1]);
f = fittype('a*exp(-x/b)+c','options',s);

[c2,gof2] = fit((0:length(soficor)-1)',mitrace(1:end)',f);
a = c2.a;
b = c2.b;
c = c2.c;
x=(0:length(soficor)-1);
fitcurve = a*exp(-x/b)+c;

results.ton.mitrace = mitrace;
results.ton.fitinfo = c2;
results.ton.fitcurve = fitcurve;
results.ton.xaxis = x;

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