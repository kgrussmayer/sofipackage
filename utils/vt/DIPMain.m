%INIT
clear all;
close all;
clc;

S.ImageNames = {
              'GreenND02_405ND01_365ND11_DU897_BV_1544',...   % Sample 1
              'GreenND02_405ND05_365ND09_DU897_BV_1561'};     % Sample 2
 

% IMAGE SETTINGS
SampleNumber = 1;                  % Image number to process [1,2]
S.WSize = 500;                     % Subsequence SOFI image size [Def. 500]
S.ImStart = 500;                   % Start of image stack (Total num. of frames must be ImStart + NumImages)
S.NumImages = 2000;                % Number of images in the stack for SOFI processing [Def. 2000,Max. 4000 (+) ]
S.PixLims = 200;                   % Maximally 200x200 central pixels [Up to about 350]
S.MaxCorrSamp = 50;                % Maximum numer of correlation samples to show [50-100 Recommended]
S.DCCoef = 0.0;                    % Added DC after filtration [Not needed]
Threads = 2;                       % Number of working threads in Parallel (Req. Par.Com.Toolb)
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
OrderSweepFIR = 50:2:80;              % FIR - Orders to sweep through
OrderSweepAVG = 100:20:300;           % FIR Moving average - Orders to sweep through
CfsIIR = linspace(0.01,0.1,10);       % List of cut-off frequencies for Filter Sweeping IIR
CfsFIR = linspace(0.05,0.1,10);       % List of cut-off frequencies for Filter Sweeping FIR
Px = 10;                              % Radius in pixels around the widefield's maximum. Filter Optimalization is calculated on this piece of data.

% BANK OF FILTERS SETTINGS
S.CFs = 0.0;        % Frequency Start [0-->1]
S.CFa = 0.2;        % Filter Bandwidth [0-->1]
S.CFo = 0.5;        % Frequency Overlap [0-->1] (Ex: 0.5 = 50% of CFa)

% SOFI AND DECONV SETTINGS
S.NumIter = 10;                  % Number of iterations used for post processing 
S.Orders = 1:4;                  % SOFI Orders to calculate
S.Enhance = 1;                   % Enhance the final SOFI Image with medfilt,wiener2 and imadjust;
S.Deconv = 'RL';                  % RL / W  / N (Richardson-Lucy / Wiener /  No Deconv


% SHOW AND CALC SETTING
S.Wiener = 1;                 % Basic preprocessing noise removal
S.SC = 0;                     % Calculate PSD, Correlation and Phase Spectrum
S.SCSpec = 1;                 % SubOption to calculate PSD
S.SCPhase = 1;                % SubOption to calculate phase Spectrum
S.SCCorr = 0;                 % SubOption to calculate Correlation

S.SaveFigures = 0;            % Save figures to a selected location
S.ShowFilt = 0;               % Show properties of the filter
S.ShowMITrace = 1;            % Show mean traces
S.ShowFilterBank = 0;         % Calculate and show bank of filters
S.ShowAcorFiltTypes = 0;      % Show Autocorrelation of different filter types
S.SweepFIRWind = 0;           % Sweep through FIR Windows
S.ShowWidefield = 0;          % Show Widefield Images of the first 2 samples

% SWEEP ALL FILTERS (DIPAllFilter function only)
S.FiltOrderFir = 100;
S.FiltOrderIir = 3;


S.ImagePath = horzcat('C:\SOFIData',filesep);          % IMAGE dir with .tif files 
S.SavePath = horzcat('D:\DIPDVD\SaveDir',filesep);     % Optional save path for all images
addpath(genpath(horzcat('D:\DIPDVD\Data',filesep)));   % Directory with SOFI routines


S.ImName = S.ImageNames{SampleNumber};
S.ImNameSamp = horzcat('Sample',num2str(SampleNumber));

%%% INITAL LOAD OF THE IMAGE %%%
[Stack,StackO] = DIPCustomImageLoad(S);
   
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



% UNCOMMENT TO NORMALIZE STACK (Produces artifacts)
[S.ImY,S.ImX,S.ImZ] = size(Stack);
% for i = 1:S.NumImages
%     Stack(:,:,i) = Stack(:,:,i)./max(max(Stack(:,:,i)));
% end
% for i = 1:S.NMax
%     StackO(:,:,i) = StackO(:,:,i)./max(max(StackO(:,:,i)));
% end


% PREPARE BANK FILTERING
if(S.ShowFilterBank)
    Olap = S.CFa*S.CFo;
    S.NF = 1 + ceil((1 - S.CFs - S.CFa)/(S.CFa - Olap));
    S.CF1 = zeros(1,S.NF);
    S.CF2 = zeros(1,S.NF);
    for i = 1:S.NF
        S.CF1(i) = S.CFs + (i-1)*S.CFa - Olap*(i-1);
        S.CF2(i) = S.CFs + i*S.CFa - Olap*(i-1);
    end
    % FIX LAST BANK FILTER
    if(S.CF2(S.NF) > 1)
        S.CF2(S.NF) = 1;
    end
end

% CORRECTION DECISION
S.CorrectionType = lower(S.CorrectionType);
switch S.CorrectionType
    case 'fir'
        S.CorrectionType = 'fir';
    case 'iir'
        S.CorrectionType = 'iir';
    case 'avg'
        S.CorrectionType = 'avg';
    otherwise
        error('Unknown Correction Type');
end

%% AWGN NOISE REMOVAL
if(S.Wiener)
    NWN = 3;
    parfor(i = 1:S.ImZ,Threads)
        Stack(:,:,i) = wiener2(squeeze(Stack(:,:,i)),[NWN,NWN]);
    end
    parfor(i = 1:S.ImZ + S.ImStart)
        [StackO(:,:,i),SNR(i)] = wiener2(squeeze(StackO(:,:,i)),[NWN,NWN]);
    end
    
      % SHOW ESTIMATED SNR
%     figure('name','Estimated Image SNR','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
%     plot(SNR,'color',S.Cl1,'linewidth',S.Lw);
%     grid on;
%     xlabel('Frame [-]','FontSize',S.FontSize);
%     ylabel('Estimated noise [\sigma^{2} ]','FontSize',S.FontSize);
        
end

if(isequal(S.ShowFilterBank,1))
    
    % DEMO ...
    
    %% TRY SYMMETRIC FILTER BANK DECOMPOSITION
    %Nfilt = 8; FiltOrder = 150;clc;
    %[ BankSignals, fig ] = DIPBankDecompositionSym( Stack, S, Nfilt, FiltOrder, Threads );

    %% TRY ASYMMETRIC FILTER BANK DECOMPOSITION
    %Nfilt = 4; FiltOrder = 150;clc;
    %[ BankSignals, fig ] = DIPBankDecompositionASym( Stack, S, Nfilt, FiltOrder, Threads );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIND FIR OPTIMUM CUT-OFF FREQUNCY AND WINDOW TYPE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(S.CorrectionType,'fir'))

if(mod(S.FiltOrder,2)== 0)
    NWin = S.FiltOrder + 1;
else
    NWin = S.FiltOrder + 2;
end

% DEFINE WINDOWS FOR FIR SWEEPING
Windows = []; Counter = 1;
Windows(Counter,:) = hamming(NWin)'; WinT{Counter} = 'Hamming'; Counter = Counter + 1;
Windows(Counter,:) = nuttallwin(NWin)'; WinT{Counter} = 'Nuttal'; Counter = Counter + 1;
Windows(Counter,:) = kaiser(NWin,0.01)'; WinT{Counter} = 'Kaiser 0.01'; Counter = Counter + 1;
Windows(Counter,:) = kaiser(NWin,0.5)'; WinT{Counter} = 'Kaiser 0.50'; Counter = Counter + 1;
Windows(Counter,:) = kaiser(NWin,0.99)'; WinT{Counter} = 'Kaiser 0.99'; Counter = Counter + 1;
Windows(Counter,:) = blackmanharris(NWin)'; WinT{Counter} = 'BlackmanHarris'; Counter = Counter + 1;
Windows(Counter,:) = gausswin(NWin)'; WinT{Counter} = 'Gauss'; Counter = Counter + 1;
Windows(Counter,:) = flattopwin(NWin)'; WinT{Counter} = 'FlatTop'; Counter = Counter + 1;

[S.CfOpt,S.WinOpt,S.FiltOrderOpt,Correl,IWin] = DIPFIRSweep( StackO, S , CfsFIR, Windows, WinT, OrderSweepFIR, Px);
disp(horzcat('FIR Optimum Cf: ',num2str(S.CfOpt)));
disp(horzcat('FIR Optimum Order: ',num2str(S.FiltOrderOpt)));
disp(horzcat('FIR Optimum Window: ',num2str(WinT{IWin})));
[N,~] = size(Correl);
Cls = jet(N);

f0 = figure('name','Find Optimum','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
for i = 1:N
    plot(Correl(i,:),'color',Cls(i,:),'linewidth',S.Lw);
    hold on;grid on;
end
xlabel('Frame [-]','FontSize',S.FontSize);
ylabel('Correlation [-]','FontSize',S.FontSize);

if(isequal(S.SaveFigures, 1))
    print(f0,horzcat(S.SavePath,S.ImName,'CFSweep.png'),'-dpng');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECK AUTOCORRELATION BASED ON USED FILTER %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(S.ShowAcorFiltTypes)
    [ Correlations,FilterTypes,FilterOrders ] = DIPAllFilter(StackO,S );
    f0 = figure('name','Different Filter Types','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    Cls = jet(length(Correlations(:,1)));
    for i = 1:length(Correlations(:,1))
        plot(Correlations(i,:),'color',Cls(i,:),'linewidth',S.Lw);
        hold on;grid on;
        Legends{i} = horzcat(num2str(FilterOrders(i)),' Order ',FilterTypes{i});
    end
    xlabel('Frame [-]','FontSize',S.FontSize);
    ylabel('Correlation [-]','FontSize',S.FontSize);
    title('Influence of filter type on the autocorrelation function','FontSize',S.FontSize);
    legend(Legends);
    
    if(isequal(S.SaveFigures, 1))
        print(f0,horzcat(S.SavePath,S.ImName,'DifferentFilterSweep.png'),'-dpng');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWEEP IIR FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(S.CorrectionType,'iir'))
    %OrderSweepIIR = 1:2;
    %CfsIIR = linspace(0.01,0.1,10);
    %S.IIRFiltTypes = {'butter','cheby','cauer'};
    [ Correlations,S.FiltOrderOpt,S.FiltType,S.CfOpt Ords, Cfss, Filts,IndOpt] = DIPIIRSweep( StackO, S, OrderSweepIIR, S.IIRFiltTypes, CfsIIR, Px);
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWEEP MOVING AVERAGE ORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(S.CorrectionType,'avg'))
    [ S.FiltOrderOpt, Correlations] = DIPAVGSweep( StackO, S ,OrderSweepAVG, Px, Threads );
    
    f0 = figure('name','Moving average','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    Cls = jet(length(OrderSweepAVG));
    for i = 1:length(OrderSweepAVG)
        plot(Correlations(i,:),'color',Cls(i,:),'Linewidth',S.Lw);
        hold on;grid on;
        Legends{i} = horzcat('Order: ',num2str(OrderSweepAVG(i)));
    end
    xlabel('Frame [-]','FontSize',S.FontSize);
    ylabel('Correlation [-]','FontSize',S.FontSize);
    legend(Legends);
    clear Legends;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWEEP FIR WINDOW TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(S.SweepFIRWind)
    S.FiltCutOff1 = 0.05;S.FiltOrderS = 100;
    [Correlations,Windows,WinT ] = DIPSweepFirWindTypes( StackO, S, S.FiltOrderS, Px, Threads);
    N = length(Correlations(:,1)); Cls = jet(N);
    
    f0 = figure('name','FIR Window Types','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    subplot(1,2,1);
    for i = 1:N
        plot(Correlations(i,:),'color',Cls(i,:),'linewidth',S.Lw);
        hold on;grid on;
        Legends{i} = horzcat(WinT{i},' Window');
    end
     xlabel('Frame [-]','FontSize',S.FontSize);
     ylabel('Correlation [-]','FontSize',S.FontSize);
     title('(a)','FontSize',S.FontSize);
     legend(Legends);
     
     S.FiltCutOff1 = 0.1;
     [Correlations,Windows,WinT ] = DIPSweepFirWindTypes( StackO, S, S.FiltOrderS, Px, Threads);
     subplot(1,2,2);
     for i = 1:N
        plot(Correlations(i,:),'color',Cls(i,:),'linewidth',S.Lw);
        hold on;grid on;
        Legends{i} = horzcat(WinT{i},' Window');
    end
     xlabel('Frame [-]','FontSize',S.FontSize);
     ylabel('Correlation [-]','FontSize',S.FontSize);
     title('(b)','FontSize',S.FontSize);
     legend(Legends);
     clear Legends;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEANTRACE & BLEACHING CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IMAGE NAME
switch S.CorrectionType
    case 'fir'
        IMName = horzcat('PIX',num2str(S.PixLims),'NUMIM',num2str(S.NumImages),'FIR',num2str(S.FiltOrderOpt),'CF',num2str(S.CfOpt),S.ImName);
    case 'iir'
        IMName = horzcat('PIX',num2str(S.PixLims),'NUMIM',num2str(S.NumImages),'IIR',S.FiltType,num2str(S.FiltOrderOpt),'CF',num2str(S.CfOpt),S.ImName);
    case 'avg'
        IMName = horzcat('PIX',num2str(S.PixLims),'NUMIM',num2str(S.NumImages),'MovingAVG',num2str(S.FiltOrderOpt),'CF',S.ImName);
end
if(S.Wiener)
    IMName = horzcat('Wiener',IMName);
end

% CUSTOM FILTER CORRECTION
switch S.CorrectionType
    case 'fir'
        [~,~,Ht,~,Nh,Fig] = DIPFilterConstruction('high','fir',S.CfOpt,[],S.FiltOrderOpt,S.FiltPassRip,S.FiltStopRip,S.ShowFilt,S.WinOpt);
    case 'iir'
        [~,~,Ht,~,Nh,Fig] = DIPFilterConstruction('high',S.FiltType,S.CfOpt,[],S.FiltOrderOpt,S.FiltPassRip,S.FiltStopRip,S.ShowFilt);
    case 'avg'
        disp(horzcat('Substract moving average,order: ',num2str(S.FiltOrderOpt)));
        af = 1; bf = ones(1,S.FiltOrderOpt) * 1/S.FiltOrderOpt;
        Ht = impz(bf,af);
        Nh = length(Ht);
end


% FASTER IMPLEMENTATION
[Ys,Xs,Zs] = size(StackO);
[Y,X,Z] = size(Stack);

if(strcmp(S.CorrectionType,'avg'))
   
    % PIXEL-INDIVIDUAL
    StackFilt = StackO;
    parfor (i = 1:Ys,Threads)
        for j = 1:Xs
            OSignal = squeeze(StackO(i,j,:))';
            MeanSig = filter(bf,af,OSignal);
            MeanSig = circshift(MeanSig,[0,-S.FiltOrderOpt/2]);
            switch S.AVGCorrectionType
                case 'S'
                    StackFilt(i,j,:) =  (OSignal- MeanSig);
                case 'I'
                    StackFilt(i,j,:) =  OSignal./MeanSig;
                case 'M'
                    StackFilt(i,j,:) =  (OSignal- MeanSig) + OSignal./MeanSig;
                otherwise
                    error('Substract AVG Correction,Please Speicify Single(S)/Inverse(I)/Multi(M)');
            end
            
        end
    end
    
    % CUT AWAY
    StackFilt = 10*StackFilt(:,:,S.ImStart-S.FiltOrderOpt/2:S.ImStart + S.NumImages-S.FiltOrderOpt/2-1);
    Stack = StackO(:,:,S.ImStart-S.FiltOrderOpt/2:S.ImStart + S.NumImages-S.FiltOrderOpt/2-1);
else
    Nfft = 2^(ceil(log2(Nh+Zs-1)));
    Nfft = 2*Nh + Zs + 1;
    Hzp = fft([Ht',zeros(1,Nfft-Nh)]);
    StackFilt = zeros(Y,X,Nfft);
    StackO = cat(3,StackO,zeros(Y,X,Nfft-Zs));
    StackFilt = fft(StackO,[],3);
    for i =1:Nfft
        StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
    end
    StackFilt = real(ifft(StackFilt,[],3));
    
    % CUT AWAY CONVOLVED DATA,RESTORE STACKO
    StackFilt = StackFilt(:,:,S.ImStart:S.ImStart + S.NumImages - 1);
    StackO = StackO(:,:,1:Zs);
end


% INVERSE EXPONENTIAL FITTING
MI = squeeze(mean(mean(Stack,1),2));
MI = MI./max(MI);


%% RUN OVERNIGHT THIS TYPE OF CORRECTION!!!
% OR TAKE A COFFEE,GO OUT FOR A WALK ...
% tic;
% [Y,X,Z] = size(Stack);
% Stack = double(Stack);
% StackExp = Stack;
% parfor(i = 1:Y,4)
%     disp(horzcat('Processing: ',num2str(i),' / ',num2str(Y)));
%     FitOpt = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0],'Upper',[Inf,Inf,Inf],'Startpoint',[1 1 1]);
%     FitType = fittype('a*exp(-x/b)+c','options',FitOpt);
%     for j = 1:X
%         Sig = squeeze(Stack(i,j,:));
%         [C2,~] = fit((0:Z-1)',Sig,FitType);
%         for k = 1:Z
%             Sig(k) = Sig(k)/(C2.a*exp(-(k-1)/C2.b)+C2.c);
%         end
%         StackExp(i,j,:) = Sig;
%     end  
% end
% tfit = toc;
% disp(horzcat('PixelWise Fit elapsed time: ',num2str(tfit),' s'));

% THIS IS NOT AS TIME CONSUMING
FitOpt = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[Inf,Inf,Inf],...
               'Startpoint',[1 1 1]);
FitType = fittype('a*exp(-x/b)+c','options',FitOpt);                    % FIT DATA
[C2,~] = fit((0:S.NumImages-1)',MI,FitType);               % FIT DATA    
disp(horzcat('Fit Coeficients a: ',num2str(C2.a),' b: ',num2str(C2.b),' c: ',num2str(C2.c)));


StackExp = Stack;
for (i = 1:S.NumImages) 
    StackExp(:,:,i) = Stack(:,:,i)./(C2.a*exp(-(i-1)/C2.b) + C2.c);
end

%%

% CALCULATE MEANTRACE AFTER CORRECTION
StackExp = StackExp./(max(StackExp(:)));
StackFilt = StackFilt./(max(StackFilt(:)));
StackFilt = StackFilt + S.DCCoef;

% SCALE DATA TO A DIFFERENT RANGE (a,b) % OPTIONAL
% [Ys,Xs,Zs] = size(StackFilt);
% a = 1;b = 2;tic;
% for i = 1:Ys
%     for j = 1:Xs
%         Sig = squeeze(StackFilt(i,j,:))';
%         Mi = min(Sig); Ma = max(Sig);
%        for k = 1:Zs
%            Sig(k) = (b-a)*(Sig(k)-Mi)/(Ma-Mi) + a;
%        end
%        StackFilt(i,j,:) = Sig;
%     end
% end
% TimeTaken = toc;
% disp(horzcat('Scaling to interval: ',num2str(a),'-',num2str(b),' Finsihed in: ',num2str(TimeTaken),' s'));
% clear Ys;
% clear Xs;
% clear Zs;

if(isequal(S.ShowMITrace,1))
    MIFilt = squeeze(mean(mean(StackFilt,1),2));
    MIFilt = MIFilt./(max(MIFilt));
    MIExp = squeeze(mean(mean(StackExp,1),2));
    MIExp = MIExp./max(MIExp);
    XAxis = (0:S.NumImages-1);
    Legends = [];
    f1 = figure('name','MEAN TRACE','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    plot(XAxis,MI,'color',S.Cl1,'linewidth',S.Lw);
    Legends{1} = 'Mean Trace - Original Image';
    hold on;grid on;
    
    plot(XAxis,MIExp,'color',S.Cl3,'linewidth',S.Lw);
    Legends{2} = 'Mean Trace after inverse fitting';
    
    plot(XAxis,MIFilt,'color',S.Cl4,'linewidth',S.Lw);
    Legends{3} = 'Mean Trace after filtration';
    
    plot(XAxis,C2.a*exp(-XAxis/C2.b)+C2.c,'color',S.Cl2,'linewidth',S.Lw);
    Legends{2} = 'Exponential Fit Curve';
    
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

[X] = DIPCalculation(Stack,StackFilt,StackExp,S);

Spec = X.Spec; Phase = X.Phase; Corela = X.Corela;
SpecFilt = X.SpecFilt; PhaseFilt = X.PhaseFilt; CorelaFilt = X.CorelaFilt;
SpecExp = X.SpecExp; PhaseExp = X.PhaseExp; CorelaExp = X.CorelaExp;


PlotSpec = 20*log10(fftshift(abs(Spec)));            % Shift and standard spectrum 20log10
PlotSpecFilt = 20*log10(fftshift(abs(SpecFilt)));    % Shift and standard spectrum 20log10
PlotSpecExp = 20*log10(fftshift(abs(SpecExp)));      % Shift and standard spectrum 20log10
SpecAx = linspace(-1,1,S.ImZ);                       % Spectrum Axis
CorelAx = 1:S.MaxCorrSamp;
clear X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOFI HR IMAGE CALCULATION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOFIO - Original | SOFI - WITHOUT DECONV | SOFIP - Deconvolved

[SOFIO,fwhms] = DIPSOFICumulants(Stack,S);
if(mean(fwhms) <1 || mean(fwhms)>5)
    warning(horzcat('FWHM Without: ',num2str(mean(fwhms))));
end
[SOFI] =  sofiLinearize(SOFIO,fwhms,[],S.Orders,S.NumIter,'N');
[SOFIP ] = sofiLinearize(SOFIO,fwhms,[],S.Orders,S.NumIter,S.Deconv);
clear SOFIO;

[SOFIFiltO,fwhms] = DIPSOFICumulants(StackFilt,S);
if(mean(fwhms) <1 || mean(fwhms)>5)
    warning(horzcat('FWHM Filter: ',num2str(mean(fwhms))));
end
[SOFIFilt] = sofiLinearize(SOFIFiltO,fwhms,[],S.Orders,S.NumIter,'N');
[SOFIFiltP] = sofiLinearize(SOFIFiltO,fwhms,[],S.Orders,S.NumIter,S.Deconv);
clear SOFIFiltO;

[SOFIExpO,fwhms] = DIPSOFICumulants(StackExp,S);
if(mean(fwhms) <1 || mean(fwhms)>5)
    warning(horzcat('FWHM ExpFit: ',num2str(mean(fwhms))));
end
[SOFIExp] = sofiLinearize(SOFIExpO,fwhms,[],S.Orders,S.NumIter,'N');
[SOFIExpP] = sofiLinearize(SOFIExpO,fwhms,[],S.Orders,S.NumIter,S.Deconv);
clear SOFIExpO;


% MULTIFRAME ORIGINAL SOFI
[SOFIMFO,fwhms] = DIPSOFICumulantsFramed(StackExp,S);
if(mean(fwhms) <1 || mean(fwhms)>5)
    warning(horzcat('FWHM ExpFitFrames: ',num2str(mean(fwhms))));
end
[SOFIMF] = sofiLinearize(SOFIMFO,mean(fwhms),[],S.Orders,S.NumIter,'N');
[SOFIMFP] = sofiLinearize(SOFIMFO,mean(fwhms),[],S.Orders,S.NumIter,S.Deconv);
clear SOFIMFO;


%% CALCULATE WITH BANK OF FILETRS
if(S.ShowFilterBank)
    [BankSofi,Correlations]  = DIPSOFIBankFilter( StackO, S );

    fa = figure('name','Filter Bank Correlations','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    Cls = hsv(S.NF);
    for i = 1:S.NF
        plot(Correlations(i,:),'color',Cls(i,:),'linewidth',S.Lw);
        hold on;grid on;
        Legends{i} = horzcat('From: ',num2str(S.CF1(i)),' To: ',num2str(S.CF2(i)));
    end
    title('Filter Bank Correlations','FontSize',S.FontSize);
    ylabel('Correlation [-]','FontSize',S.FontSize);
    xlabel('Frame [-]','FontSize',S.FontSize);
    legend(Legends);
end


%% SHOW WIDEFIELD OF THE TWO SAMPLES
if(S.ShowWidefield)
    S.ImName = S.ImageNames{1};
    StackW1 = DIPCustomImageLoad(S);
    S.ImName = S.ImageNames{2};
    StackW2 = DIPCustomImageLoad(S);
    Widefield1 = squeeze(mean(StackW1,3));
    Widefield2 = squeeze(mean(StackW2,3));
    IMS{1} = Widefield1;
    IMS{2} = Widefield2;
    Titles = {'(a)','(b)'};
    
    figure('name','Widefield Image','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    
    Nr = 1; % NumRows
    Nc = 2; % NumColumns
    Spacing = 0.03;
    Wis = (1-(Nc+1)*Spacing)/Nc ;
    Hes = (1-(Nr+1)*Spacing)/Nr ;

    % PLOT IMAGES
    for i = 1:Nr
        for j = 1:Nc
            Ind = j + (i-1)*Nc;

            W = (j-1)*Wis + (j)*Spacing ;
            H = (i-1)*Hes + (i)*Spacing ;
            hx = axes('position',[W,H,Wis,Hes]);h(Ind) = hx;

            IM = IMS{Ind};
            imshow(IM,[min(IM(:)),max(IM(:))]);
            title(Titles{Ind},'FontSize',S.FontSize);

        end
    end
    colormap morgenstemning;
    clear StackW1;
    clear StackW2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHOW SOFI IMAGES %%%%%%%
% SET DIFFERENT WINDOW SIZE 
if(S.ShowFilterBank)
    S.Wi = 1000;
    S.He = round(S.He*2/3);
else
    S.He = 1000;
    S.Wi = 1000;
end

NMed = 5;

for ind = 2:length(S.Orders)

    switch S.CorrectionType
        case 'fir'
            IMNameS = horzcat(S.ImNameSamp,'SOFIOrd',num2str(ind),'Pix',num2str(S.PixLims),'Frames',num2str(S.NumImages),'FIR',num2str(S.FiltOrderOpt),'CF',num2str(S.CfOpt));
        case 'iir'
            IMNameS = horzcat(S.ImNameSamp,'SOFIOrd',num2str(ind),'Pix',num2str(S.PixLims),'Frames',num2str(S.NumImages),'IIR',S.FiltType,num2str(S.FiltOrderOpt),'CF',num2str(S.CfOpt));
        case 'avg'
            IMNameS = horzcat(S.ImNameSamp,'SOFIOrd',num2str(ind),'Pix',num2str(S.PixLims),'Frames',num2str(S.NumImages),'MovingAVG-',S.AVGCorrectionType,num2str(S.FiltOrderOpt));
    end

    if(S.Wiener)
        IMNameS = horzcat(IMNameS,'WienerDenoise');
    end
    

    %%% WITHOUT DECONVOLUTION
    fs(ind) = figure('name',horzcat('SOFI Order: ',num2str(ind)),'position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);

    %%% ONLY EXP-FIT AND DIGITAL FILTER
    if(S.ShowFilterBank)
        subplot(1,3,1);
        IM = mean(SOFIMF{ind},3); % SOFIMFP for deconv result
        if(S.Enhance)
            IM = (medfilt2((( IM )),[ NMed, NMed]));
            IM = wiener2(IM,[NMed,NMed]);
        end
        IM = imadjust(IM,stretchlim(IM),[]);
        imshow(IM,[min(IM(:)),max(IM(:))]);
        title('(a)','FontSize',S.FontSize);


        subplot(1,3,2);
        IM = SOFIFilt{ind}; % SOFIFiltP for deconv result
        if(S.Enhance)
            IM = (medfilt2((( IM )),[ NMed, NMed]));
            IM = wiener2(IM,[NMed,NMed]);
        end
        IM = imadjust(IM,stretchlim(IM),[]);
        imshow(IM,[min(IM(:)),max(IM(:))]);
        title('(b)','FontSize',S.FontSize);

        subplot(1,3,3);
        IM = BankSofi{ind};
        if(S.Enhance)
            IM = (medfilt2((( IM )),[ NMed, NMed]));
            IM = wiener2(IM,[NMed,NMed]);
        end
        IM = imadjust(IM,stretchlim(IM),[]);
        imshow(IM,[min(IM(:)),max(IM(:))]);
        title('(c)','FontSize',S.FontSize);
        colormap morgenstemning;
    else
        
        % a ... Without Correction
        % b ... Inverse Fit
        % c ... Multiframe Inverse fit
        % d ... Filter
        
        % PREPARE IMAGES
        IMS{1} = SOFI{ind};     % SOFIP for deconv result
        IMS{2} = SOFIExp{ind};  % SOFIExpP for deconv result
        IMS{3} = SOFIMF{ind};   % SOFIMFP for deconv result
        IMS{4} = SOFIFilt{ind}; % SOFIFiltP for deconv result
        for i = 1:4
            if(S.Enhance)
                IM = IMS{i};
                [~,~,Z] = size(IM);
                if(Z>1)
                    IM = mean(IM,3);
                end
                
                IM = (medfilt2((( IM )),[ NMed, NMed]));
                IM = wiener2(IM,[NMed,NMed]);
            end
            IMS{i} = imadjust(IM,stretchlim(IM),[]);
            
        end
        Titles = {'(a)','(b)','(c)','(d)'};
        
        
        
        
        Nr = 2; % NumRows
        Nc = 2; % NumColumns
        Spacing = 0.03;
        Wis = (1-(Nc+1)*Spacing)/Nc ;
        Hes = (1-(Nr+1)*Spacing)/Nr ;
        
        % PLOT IMAGES
        for i = 1:Nr
            for j = 1:Nc
                Ind = j + (i-1)*Nc;

                W = (j-1)*Wis + (j)*Spacing ;
                H = (Nc-i)*Hes + (Nc - i+1)*Spacing ;
                hx = axes('position',[W,H,Wis,Hes]);h(Ind) = hx;

                IM = IMS{Ind};
                imshow(IM,[min(IM(:)),max(IM(:))]);
                title(Titles{Ind},'FontSize',S.FontSize);

            end
        end
        colormap morgenstemning;

    end
    
    % DECONV - NOT DECONV COMPARISON
    fsd(ind) = figure('name','Not Deconv,Deconv','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He/2]);
    IMS{1} = SOFIFilt{ind};
    IMS{2} = SOFIFiltP{ind};
    
    for i = 1:2
        if(S.Enhance)
            IM = IMS{i};
            [~,~,Z] = size(IM);
            if(Z>1)
                IM = mean(IM,3);
            end

            IM = (medfilt2((( IM )),[ NMed, NMed]));
            IM = wiener2(IM,[NMed,NMed]);
        end
        IMS{i} = imadjust(IM,stretchlim(IM),[]);

    end
    Titles = {'(a)','(b)'};

    Nr = 1; % NumRows
    Nc = 2; % NumColumns
    Spacing = 0.03;
    Wis = (1-(Nc+1)*Spacing)/Nc ;
    Hes = (1-(Nr+1)*Spacing)/Nr ;

    % PLOT IMAGES
    for i = 1:Nr
        for j = 1:Nc
            Ind = j + (i-1)*Nc;

            W = (j-1)*Wis + (j)*Spacing ;
            H = (i-1)*Hes + (i)*Spacing ;
            hx = axes('position',[W,H,Wis,Hes]);h(Ind) = hx;

            IM = IMS{Ind};
            imshow(IM,[min(IM(:)),max(IM(:))]);
            title(Titles{Ind},'FontSize',S.FontSize);

        end
    end
    colormap morgenstemning;





    if(S.SaveFigures)
        saveas(fsd(ind),horzcat(S.SavePath,'DECONV',IMNameS,'.fig'),'fig');
        saveas(fs(ind),horzcat(S.SavePath,IMNameS,'.fig'),'fig');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHOW SPECTRUM OF TIME SIGNALS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isequal(S.SC,1) && ( isequal(S.SCSpec,1) || isequal(S.SCPhase,1)) )
    f2 = figure('name','Spectrum','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
    if(S.SCSpec)
        subplot(2,1,1);
        plot(SpecAx,PlotSpec,'color',S.Cl1,'linewidth',S.Lw);
        hold on;grid on;
        plot(SpecAx,PlotSpecExp,'color',S.Cl2,'linewidth',S.Lw);
        plot(SpecAx,PlotSpecFilt,'color',S.Cl3,'linewidth',S.Lw);
        legend('Without Correction','Inv. exp. fit ','Filtration');
        ylabel('20log_{10}(|S(\omega)|)','FontSize',S.FontSize);
        xlabel('Normalized Frequency x\pi rad/sample','FontSize',S.FontSize);
        %title('Average Spectrum','FontSize',S.FontSize);
    end
    
    
    
    if(S.SCPhase)
        subplot(2,1,2);
        plot(SpecAx,Phase,'color',S.Cl1,'linewidth',S.Lw);
        hold on;grid on;
        plot(SpecAx,PhaseExp,'color',S.Cl2,'linewidth',S.Lw);
        plot(SpecAx,PhaseFilt,'color',S.Cl3,'linewidth',S.Lw);
        legend('Without correction','Inv. exp. fit ','Filtration');
        xlabel('Normalized frequency x \pi rad/sample','FontSize',S.FontSize);
        ylabel('Phase [°deg]','FontSize',S.FontSize);
    end

    if(S.SaveFigures)
        print(f2,horzcat(S.SavePath,'Spectrum',IMName,'.png'),'-dpng');
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHOW NORMALIZED CORRELATION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isequal(S.SC,1) && isequal(S.SCCorr,1))
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
    
    if(S.SaveFigures)
        print(f3,horzcat(S.SavePath,'Correlation',IMName,'.png'),'-dpng');
    end
end





