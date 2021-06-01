function [ CFOptimal,WinOptimal,OrdOptimal,Correlations, Wind ] = DIPFIRSweep( Stack, S , Cfs, Windows, WinT, Orders, Px)
% Stack ...     Input Image
% S ...         Input Parameters
% Cfs ...       Vector of Selected cut-off frequencies
% Windows ...   Windows for optimum FIR Design
% WinT ...      Window names
% Orders ...    Orders to sweep through
% Px ...        PixelSize radius for quick acquisition
%
%
% OUT:
% CFOptimal ... optimum CutOff frequency for a given Hug-Pass Filter

CFOptimal = 0.5;
Nf = length(Cfs);
Nw = length(WinT);
No = length(Orders);
MaxSamp = S.MaxCorrSamp;

% FIND MAXIMUM
Widefield = sum(Stack,3);
[M,I] = max(Widefield);
[MM,II] = max(M);
Yi = I(II);
Xi = II;

% FIX EDGES PROBLEMS
[Y,X] = size(Widefield);
if(Xi-Px <=0 )
    Xi = Xi + Px;
end
if(Xi + Px > X)
    Xi = Xi - Px;
end
if(Yi - Px <= 0)
    Yi = Yi + Px;
end
if(Yi + Px >= Y)
    Yi = Yi - Px;
end

% SEARCH FOR OPTIMUM IN THIS ARRAY INSTEAD
StackC = Stack(Yi-Px:Yi+Px , Xi-Px:Xi+Px , :);


% COMPUTING PART
Correlations = zeros(Nf*Nw,S.MaxCorrSamp);

% SIMPLIFY INDEXING
TempCfs =   zeros(1,Nf * Nw);                  % Cut-Off Frequencies
TempWins =  zeros(Nf*Nw,length(Windows(1,:))); % Windows (Samples)
TempWinsT = cell(1,Nf*Nw);                     % Window Names

Counter = 1;
for iNf = 1:Nf
    for iNw = 1:Nw
        CorelaFilt = [];
        StackFilt = [];
        
        % LINEAR INDEXING
        Ind = iNw + (iNf-1)*Nw; 
        disp(horzcat('FIR Opt Processing: ',num2str(Counter),' / ',num2str(Nf*Nw)));
        
        % CONSTRUCT FIR FILTER
        Cf = Cfs(iNf);TempCfs(Ind) = Cf;
        Win = Windows(iNw,:);TempWins(Ind,:) = Win;
        TempWinsT(Ind) = WinT(iNw);
      
        [~,~,Ht,~,Nh,~] = DIPFilterConstruction('high','fir',Cf, 0 ,S.FiltOrder,S.FiltPassRip,S.FiltStopRip,0,Win);

        Ns = length(StackC(1,1,:)); 
        Nfft = 2^(ceil(log2(Nh+Ns-1)));
        Hzp = fft([Ht',zeros(1,Nfft-Nh)]);
        [Y,X,Z] = size(StackC);

        % FILTER STACK WITH A SELECTED FILTER
        [Yt,Xt,Zt] = size(StackC);
        StackC = cat(3,StackC,zeros(Yt,Xt,Nfft-Zt));
        StackFilt = fft(StackC,[],3);
        for  i = 1:Nfft
            StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
        end
        StackFilt = real(ifft(StackFilt,[],3));
        StackC = StackC(:,:,1:Zt);

        % CUT AFTER FILTRATION
        StackFilt = StackFilt(:,:,S.ImStart:S.ImStart + S.NumImages - 1);
        [Y,X,Z] = size(StackFilt);
        
        % COMPUTE CORRELATION
        CorelaFilt = cat(3,StackFilt,zeros(Y,X,Z));
        CorelaFilt = abs(ifft(abs(fft(CorelaFilt,[],3)).^2,[],3));
        CorelaFilt = squeeze(mean(mean(CorelaFilt,1),2));
        CorelaFilt = CorelaFilt(1:S.MaxCorrSamp);
        CorelaFilt = CorelaFilt./max(CorelaFilt);
        Correlations(Ind,:) = CorelaFilt;
        Counter = Counter + 1;
    end
    
end


%% FIND OPTIMAL FILTER BASED ON WINDOW AND CF
Differences = zeros(1,Nf*Nw);
for i = 1:Nf*Nw    
    
    % FIND LOCAL MINIMUM
    Lmin = MaxSamp; % Index Init
    for j = 2:MaxSamp
        if(Correlations(i,j) >= Correlations(i,j-1))
            Lmin = j-1;
            break;
        end
    end
    Mil = Correlations(i,Lmin);
    [Mal,Lmax] = max(Correlations(i,Lmin:end));
    PeakCoef = Mal - Mil;
    
    if(isequal(Lmin,Lmax))
        PeakCoef = 1;
    end
    
    % WEIGTENING BY MEAN VALUE
    Differences(i) = (10*PeakCoef + 1) * mean(Correlations(i,:));
end


% F1 W1 W2 W3 ... F2 W1 W2 W3 ... F3 W1 W2 W3 ...
[~,I] = min(Differences);
CFOptimal = TempCfs(I);
WinOptimal = TempWins(I,:);
Wind = mod(I,Nw)+1;


%% SWEEP FILTER ORDER
CorrelationsOrig = Correlations;
clear Correlations;
for k = 1:No
    
    Ord = Orders(k);
    if(mod(Ord,2)== 0)
        NWin = Ord + 1;
    else
        NWin = Ord + 2;
    end
    
    switch Wind
        case 1
            CurrentWindow = hamming(NWin)';
        case 2
            CurrentWindow = nuttallwin(NWin)';
        case 3
            CurrentWindow = kaiser(NWin,0.01)';
        case 4
            CurrentWindow = kaiser(NWin,0.5)';
        case 5
            CurrentWindow = kaiser(NWin,0.99)';
        case 6
            CurrentWindow = blackmanharris(NWin)';
        case 7
            CurrentWindow = gausswin(NWin)';
        case 8
            CurrentWindow = flattopwin(NWin)';
        otherwise
            error('Unknown Window Function');
    end
     [~,~,Ht,~,Nh,~] = DIPFilterConstruction('high','fir',CFOptimal, 0, Ord,[],[],[],CurrentWindow);
     
     Ns = length(StackC(1,1,:)); 
     Nfft = 2^(ceil(log2(Nh+Ns-1)));
     Hzp = fft([Ht',zeros(1,Nfft-Nh)]);
     [Y,X,Z] = size(StackC);

     % FILTER STACK WITH A SELECTED FILTER
     [Yt,Xt,Zt] = size(StackC);
     StackC = cat(3,StackC,zeros(Yt,Xt,Nfft-Zt));
     StackFilt = fft(StackC,[],3);
     for  i = 1:Nfft
         StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
     end
     
     % CUT AFTER FILTRATION
     StackFilt = real(ifft(StackFilt,[],3));
     StackC = StackC(:,:,1:Zt);
     
     % COMPUTE CORRELATION
     CorelaFilt = cat(3,StackFilt,zeros(Y,X,Z));
     CorelaFilt = abs(ifft(abs(fft(CorelaFilt,[],3)).^2,[],3));
     CorelaFilt = squeeze(mean(mean(CorelaFilt,1),2));
     CorelaFilt = CorelaFilt(1:S.MaxCorrSamp);
     CorelaFilt = CorelaFilt./max(CorelaFilt);
     Correlations(k,:) = CorelaFilt;
end

%% SELECT THE BEST FILTER ORDER
Differences = zeros(1,No);
for i = 1:No   
    
    % FIND LOCAL MINIMUM
    Lmin = MaxSamp; % Index Init
    for j = 2:MaxSamp
        if(Correlations(i,j) >= Correlations(i,j-1))
            Lmin = j-1;
            break;
        end
    end
    Mil = Correlations(i,Lmin);
    [Mal,Lmax] = max(Correlations(i,Lmin:end));
    PeakCoef = Mal - Mil;
    
    if(isequal(Lmin,Lmax))
        PeakCoef = 1;
    end
    
     disp(horzcat('Peak: ',num2str(10*PeakCoef),' Order: ',num2str(Orders(i))));
    % WEIGTENING BY MEAN VALUE
    Differences(i) = (10*PeakCoef + 1) * mean(Correlations(i,:));
end

[~,I] = min(Differences);
OrdOptimal = Orders(I);

%% FINALLY INCREASE THE FILTER ORDER
Increased = 0;
while 1
    % INCREASE ALWAYS BY 2
    if(mod(OrdOptimal + 2,2)== 0)
        NWin = OrdOptimal + 1 + 2;
    else
        NWin = OrdOptimal + 2 + 2;
    end
    if(mod(OrdOptimal,2)== 0)
        NWinO = OrdOptimal + 1;
    else
        NWinO = OrdOptimal + 2;
    end
    
    % SWITCH WINDOW
    switch(Wind)
        case 1
            WinTest = hamming(NWin)';
            WinOpt = hamming(NWinO)';
        case 2
             WinTest = nuttallwin(NWin)';
             WinOpt = nuttallwin(NWinO)';
        case 3
             WinTest = kaiser(NWin,0.01)';
             WinOpt = kaiser(NWinO,0.01)';
        case 4
             WinTest = kaiser(NWin,0.5)';
             WinOpt = kaiser(NWinO,0.5)';
        case 5
             WinTest = kaiser(NWin,0.99)';
             WinOpt = kaiser(NWinO,0.99)';
        case 6
             WinTest = blackmanharris(NWin)';
             WinOpt = blackmanharris(NWinO)';
        case 7
             WinTest = gausswin(NWin)';
             WinOpt = gausswin(NWinO)';
        case 8
             WinTest = flattopwin(NWin)';
             WinOpt = flattopwin(NWinO)';
    end
    
    % NEW FILT ORDER
    a = 1;b = fir1(OrdOptimal + 2,CFOptimal,'high',WinTest);
    [hn,~] = freqz(b,a);hn = 20*log10(hn');
    
    % CURRENT FILT ORDER
    a = 1;b = fir1(OrdOptimal,CFOptimal,'high',WinOpt);
    [hc,~] = freqz(b,a); hc = 20*log10(hc');
    
    if(hn(1) < hc(1))
        OrdOptimal = OrdOptimal +2;
        Increased = 1;
    else
        % SUFFICIENT
        break;
    end
    
end
if(Increased)
    %disp(horzcat('Filter order increased to: ',num2str(OrdOptimal)));
end

% UPDATE WINDOW FILTER LENGTH
if(mod(OrdOptimal,2)== 0)
    NWin = OrdOptimal + 1;
else
    NWin = OrdOptimal + 2;
end
switch(Wind)
    case 1
        WinOptimal = hamming(NWin)';
    case 2
         WinOptimal = nuttallwin(NWin)';
    case 3
         WinOptimal = kaiser(NWin,0.01)';
    case 4
         WinOptimal = kaiser(NWin,0.5)';
    case 5
         WinOptimal = kaiser(NWin,0.99)';
    case 6
         WinOptimal = blackmanharris(NWin)';
    case 7
         WinOptimal = gausswin(NWin)';
    case 8
         WinOptimal = flattopwin(NWin)';
end

% INCLUDE ALL RESULTS
Correlations = [CorrelationsOrig;Correlations;];


end

