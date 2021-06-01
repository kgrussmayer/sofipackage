function StackFilt = bleachcor_iir(stack,settings)
%  wraper function for various bleaching correction algorithms 

% IMAGE SETTINGS
% SampleNumber = 1;                  % Image number to process [1,2]
% S.WSize = 500;                     % Subsequence SOFI image size [Def. 500]
S.ImStart = 1;                     % Start of image stack (Total num. of frames must be ImStart + NumImages)
S.NumImages = size(stack,3);                % Number of images in the stack for SOFI processing [Def. 2000,Max. 4000 (+) ]
S.PixLims = 200;                   % Maximally 200x200 central pixels [Up to about 350]
S.MaxCorrSamp = 50;                % Maximum numer of correlation samples to show [50-100 Recommended]
S.DCCoef = 0.0;                    % Added DC after filtration [Not needed]
% Threads = 2;                       % Number of working threads in Parallel (Req. Par.Com.Toolb)
S.NMax = S.NumImages + S.ImStart;  % Maximum number of loaded Image Frames

% FILTER SETTING
S.CorrectionType = 'AVG';      % FIR / IIR / AVG - MAIN CORRECTION TYPE
S.AVGCorrectionType = 'S';     % M(Multi) - Subtract average + Inv. Fit / S(Single) Substract Mean / I(Inverse Fit)
S.FiltPassRip = 3;             % Allowed pass ripple in the pass-band (Only IIR - Cheby / Cauer)
S.FiltStopRip = 80;            % Allowed stop ripple in the stop-band (Only IIR - InvCheby / Cauer)
S.FiltCutOff1 = 0.07;          % Normalized cut-off frequency for filter design [Approx,will be replaced by optimal]
S.FiltOrder = 100;             % Prefered Filter order (FIR ~ 100 | IIR max ~ 20)
S.IIRFiltTypes = {'butter'};   %,'cheby','cauer','invcheby'}; % Butter works best

OrderSweepIIR = 1;                    % IIR - Orders to sweep through
% OrderSweepFIR = 50:2:80;              % FIR - Orders to sweep through
% OrderSweepAVG = 100:20:300;           % FIR Moving average - Orders to sweep through
CfsIIR = linspace(0.01,0.1,10);       % List of cut-off frequencies for Filter Sweeping IIR
% CfsFIR = linspace(0.05,0.1,10);       % List of cut-off frequencies for Filter Sweeping FIR
Px = 10;                              % Radius in pixels around the widefield's maximum. Filter Optimalization is calculated on this piece of data.

% BANK OF FILTERS SETTINGS
% S.CFs = 0.0;        % Frequency Start [0-->1]
% S.CFa = 0.2;        % Filter Bandwidth [0-->1]
% S.CFo = 0.5;        % Frequency Overlap [0-->1] (Ex: 0.5 = 50% of CFa)

% SOFI AND DECONV SETTINGS
% S.NumIter = 10;                  % Number of iterations used for post processing 
% S.Orders = 1:4;                  % SOFI Orders to calculate
% S.Enhance = 1;                   % Enhance the final SOFI Image with medfilt,wiener2 and imadjust;
% S.Deconv = 'RL';                  % RL / W  / N (Richardson-Lucy / Wiener /  No Deconv

% SHOW AND CALC SETTING
S.Wiener = 1;                 % Basic preprocessing noise removal
S.SC = 0;                     % Calculate PSD, Correlation and Phase Spectrum
S.SCSpec = 1;                 % SubOption to calculate PSD
S.SCPhase = 1;                % SubOption to calculate phase Spectrum
S.SCCorr = 1;                 % SubOption to calculate Correlation

S.SaveFigures = 0;            % Save figures to a selected location
S.ShowFilt = 0;               % Show properties of the filter
S.ShowMITrace = 1;            % Show mean traces
S.ShowFilterBank = 0;         % Calculate and show bank of filters
S.ShowAcorFiltTypes = 0;      % Show Autocorrelation of different filter types
S.SweepFIRWind = 0;           % Sweep through FIR Windows
S.ShowWidefield = 0;          % Show Widefield Images of the first 2 samples

% SWEEP ALL FILTERS (DIPAllFilter function only)
S.FiltOrderFir = 100;
S.FiltOrderIir = 3;                             % Radius in pixels around the widefield's maximum. Filter Optimalization is calculated on this piece of data.
 
% DEFINE COLORS
S.Cl1 = [127,0,204]/255;     % Violet
S.Cl2 = [0,128,255]/255;     % Blue
S.Cl3 = [204,0,102]/255;     % Red
S.Cl4 = [64,64,64]/255;      % Dark grey
S.Lw = 2;                    % Linewidth for plotting
S.FontSize = 13;             % FontSize for Labels and Titles

% GET DISPLAY RESOLUTION
S.Displ = groot();
S.Width = S.Displ.ScreenSize(3);
S.Height = S.Displ.ScreenSize(4);
S.Wi = 1080; % Default width of figures
S.He = 720;  % Default height of figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWEEP IIR FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OrderSweepIIR = 1:2;
%CfsIIR = linspace(0.01,0.1,10);
%S.IIRFiltTypes = {'butter','cheby','cauer'};
[ Correlations,S.FiltOrderOpt,S.FiltType,S.CfOpt, Ords, Cfss, Filts,IndOpt] = DIPIIRSweep( stack, S, OrderSweepIIR, S.IIRFiltTypes, CfsIIR, Px);
disp(horzcat('Optimum IIR Filter: ',S.FiltType,' Order: ',num2str(S.FiltOrderOpt),' CF: ',num2str(S.CfOpt)));

N = length(Correlations(:,1));
Cls = jet(N);Legends = [];

f0 = figure('name','IRR Sweep','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
for i = 1:N

    plot(Correlations(i,:),'color',Cls(i,:),'linewidth',S.Lw);
    hold on; grid on;
    Legends{i} = horzcat('Type: ',Filts{i},' Order: ',num2str(Ords(i)),' \omega_{c}: ',num2str(Cfss(i)));

end
plot(Correlations(IndOpt,:),'color',[32,32,32]/255,'linewidth',3);
xlabel('Frame [-]','FontSize',S.FontSize);
ylabel('Correlation [-]','FontSize',S.FontSize);
legend(Legends,'Optimal'); %Legends /  TOO MANY
clear Legends;

[~,~,Ht,~,Nh,Fig] = DIPFilterConstruction('high',S.FiltType,S.CfOpt,[],S.FiltOrderOpt,S.FiltPassRip,S.FiltStopRip,S.ShowFilt);

[Y,X,Zs] = size(stack);
%     [Y,X,Z] = size(Stack);

 Nfft = 2^(ceil(log2(Nh+Zs-1)));
Nfft = 2*Nh + Zs + 1;
Hzp = fft([Ht',zeros(1,Nfft-Nh)]);
StackFilt = zeros(Y,X,Nfft);
stack = cat(3,stack,zeros(Y,X,Nfft-Zs));
StackFilt = fft(stack,[],3);
for i =1:Nfft
    StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
end
StackFilt = real(ifft(StackFilt,[],3));

% CUT AWAY CONVOLVED DATA,RESTORE STACK
StackFilt = StackFilt(:,:,S.ImStart:S.ImStart + S.NumImages - 1);
stack = stack(:,:,1:Zs);

%%
% CALCULATE MEANTRACE AFTER CORRECTION
% StackExp = StackExp./(max(StackExp(:)));
StackFilt = StackFilt./(max(StackFilt(:)));
StackFilt = StackFilt + S.DCCoef;
mitrace=squeeze(mean(mean(stack,1),2));

if(isequal(S.ShowMITrace,1))
    MIFilt = squeeze(mean(mean(StackFilt,1),2));
    MIFilt = MIFilt./(max(MIFilt));
%     MIExp = squeeze(mean(mean(StackExp,1),2));
%     MIExp = MIExp./max(MIExp);
    XAxis = (0:S.NumImages-1);
    Legends = [];
    f1 = figure('name','MEAN TRACE','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    plot(XAxis,mitrace,'color',S.Cl1,'linewidth',S.Lw);
    Legends{1} = 'Mean Trace - Original Image';
    hold on;grid on;
    
%     plot(XAxis,MIExp,'color',S.Cl3,'linewidth',S.Lw);
%     Legends{2} = 'Mean Trace after inverse fitting';
    
    plot(XAxis,MIFilt,'color',S.Cl4,'linewidth',S.Lw);
    Legends{2} = 'Mean Trace after filtration';
    
%     plot(XAxis,C2.a*exp(-XAxis/C2.b)+C2.c,'color',S.Cl2,'linewidth',S.Lw);
%     Legends{2} = 'Exponential Fit Curve';
    
    ylabel('Intensity [-]','FontSize',S.FontSize);
    xlabel('Image Frame [-]','FontSize',S.FontSize);
    legend(Legends);
    clear Legends;

    if(S.SaveFigures)
       print(f1,horzcat(S.SavePath,'MeanTrace',IMName,'.png'),'-dpng');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE CORRELATION & SPECTRUM %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X] = DIPCalculation(stack(:,:,50:end),StackFilt(:,:,50:end),S);

Spec = X.Spec; Phase = X.Phase; Corela = X.Corela;
SpecFilt = X.SpecFilt; PhaseFilt = X.PhaseFilt; CorelaFilt = X.CorelaFilt;
SpecExp = X.SpecExp; PhaseExp = X.PhaseExp; CorelaExp = X.CorelaExp;
CorelAx = 1:S.MaxCorrSamp;
clear X;

Corela = Corela./max(Corela(:));
%     CorelaExp = CorelaExp./max(CorelaExp(:));
CorelaFilt = CorelaFilt./max(CorelaFilt(:));

f3 = figure('name','Correlation','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
plot(CorelAx,Corela,'color',S.Cl1,'linewidth',S.Lw);
hold on;grid on;
%     plot(CorelAx,CorelaExp,'color',S.Cl2,'linewidth',S.Lw);
plot(CorelAx,CorelaFilt,'color',S.Cl3,'linewidth',S.Lw);
legend('Without Bleaching Correction','Basic Bleaching Correction','Advanced Bleaching Correction');
xlabel('\tau','FontSize',S.FontSize);
ylabel('R(\tau)','FontSize',S.FontSize);
