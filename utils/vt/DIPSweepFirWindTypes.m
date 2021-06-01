function [ Correlations,Windows,WinT ] = DIPSweepFirWindTypes( StackO, S, Ord, Px, Threads)
% TESTING FUNCTION - SHOWING THE IMPACT OF CHOOSING SEVERAL
% DIFFERENT WINDOW FUNCTIONS FOR FIR FILTER DESIGN

MaxSamp = S.MaxCorrSamp;
Cf = S.FiltCutOff1;
ImStart = S.ImStart;
NumImg = S.NumImages;
Ns = length(StackO(1,1,:)); 
Correlations = [];


% FIND MAXIMUM
Widefield = sum(StackO,3);
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
StackC = StackO(Yi-Px:Yi+Px , Xi-Px:Xi+Px , :);

if(mod(Ord,2)== 0)
    NWin = Ord+1;
else
    NWin = Ord+2;
end
C = 1;
Windows(C,:) = hamming(NWin)'; WinT{C} = 'Hamming'; C = C+1;
Windows(C,:) = nuttallwin(NWin)'; WinT{C} = 'Nuttal';C = C+1;
Windows(C,:) = kaiser(NWin,0.25)'; WinT{C} = 'Kaiser 0.25';C = C+1;
Windows(C,:) = kaiser(NWin,0.5)'; WinT{C} = 'Kaiser 0.50';C = C+1;
Windows(C,:) = kaiser(NWin,0.75)'; WinT{C} = 'Kaiser 0.75';C = C+1;
Windows(C,:) = blackmanharris(NWin)'; WinT{C} = 'BlackmanHarris';C = C+1;
Windows(C,:) = hann(NWin)'; WinT{C} = 'Hann';C = C+1;
Windows(C,:) = gausswin(NWin)'; WinT{C} = 'Gauss';C = C+1;
Windows(C,:) = taylorwin(NWin)'; WinT{C} = 'Taylor';C = C+1;
Windows(C,:) = flattopwin(NWin)'; WinT{C} = 'FlatTop';C = C+1;
Windows(C,:) = bartlett(NWin)'; WinT{C} = 'Triangular';C = C+1;
Windows(C,:) = rectwin(NWin)'; WinT{C} = 'Rectangular';C = C+1;
Windows(C,:) = tukeywin(NWin,0.75)'; WinT{C} = 'Tukey';C = C+1;
N = C-1;

for k = 1:N
    
    disp(horzcat('Processing: ',num2str(k),' / ',num2str(N)));
    %% FILTERING 
    [~,~,Ht,~,Nh,~] = DIPFilterConstruction('high','fir',Cf,0,Ord,0,0,0,Windows(k,:));
    Nfft = 2^(ceil(log2(Nh+Ns-1)));
    Hzp = fft([Ht',zeros(1,Nfft-Nh)]);
    [Y,X,Z] = size(StackC);
    
    StackFilt = StackC;
    StackFilt = cat(3,StackFilt,zeros(Y,X,Nfft-Ns));
    StackFilt = fft(StackFilt,[],3);
    for i = 1:Nfft
        StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
    end
    StackFilt = real(ifft(StackFilt,[],3));
    StackFilt = StackFilt(:,:,S.ImStart:S.ImStart + S.NumImages - 1);

    % PROCESS AUTOCORRELATION
    StackFilt = abs(fft(cat(3,StackFilt,zeros(Y,X,NumImg)),[],3));
    StackFilt = squeeze(mean(mean(abs(ifft(StackFilt.^2,[],3)),1)));
    Correlations(k,:) = StackFilt./max(StackFilt);
    
end

% CUT CORRELATIONS
Correlations = Correlations(:,1:MaxSamp);

end

