function [ Correlations,OrderOpt,FiltTypeOpt,CfOpt,TempOrders,TempCfs,TempFilts,IndexOpt ] = DIPIIRSweep( StackO, S, Orders, FiltTypes, Cfs, Px)
% FIND THE OPTIMAL FILTER FOR FILTRATION

PassRip = S.FiltPassRip;
StopRip = S.FiltStopRip;
MaxSamp = S.MaxCorrSamp;
ImStart = S.ImStart;
NumImg = S.NumImages;

% INIT OUTPUT VARIABLES
Correlations = [];
OrderOpt = 1;
FiltTypeOpt = 'butter';
CfOpt = 0.07;

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

No = length(Orders);
Nf = length(FiltTypes);
Nc = length(Cfs);
[Y,X,Z] = size(StackC);

% FOR EASIER SEARCHING
TempOrders = zeros(1,Nf*Nc*No);
TempCfs = zeros(1,Nf*Nc*No);
TempFilts = [];
IndexOpt = 1;
Counter = 1;

for ko = 1:No % Order sweeep
    for kf = 1: Nf % filter sweep
        for kc = 1:Nc % Cut off frequency sweep
            
            % LOAD CURRENT SETTINGS
            Ord = Orders(ko);
            Cf = Cfs(kc);
            FiltType = FiltTypes{kf};
            
            % Linear Indexing
            Ind = kc + (kf-1)*Nc + (ko-1) * Nf * Nc;
            TempOrders(Ind) = Ord;
            TempCfs(Ind) = Cf;
            TempFilts{Ind} = FiltType;
            
            % FILTER CONSTRUCTION
            [ ~,~,Ht,~,Nh,~ ] = DIPFilterConstruction('high',FiltType,Cf,[],Ord,PassRip,StopRip,[] );
            Nfft = 2^(ceil(log2(Nh + Z - 1)));
            Hzp = fft([Ht',zeros(1,Nfft - Nh)]);

            % FILTER STACK
            [Yt,Xt,Zt] = size(StackC);
            StackC = cat(3,StackC,zeros(Yt,Xt,Nfft-Zt));
            StackFilt = fft(StackC,[],3);
            for i = 1:Nfft
                StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
            end
            StackFilt = real(ifft(StackFilt,[],3));
            StackC = StackC(:,:,1:Zt);
            StackFilt = StackFilt(:,:,ImStart:ImStart + NumImg - 1);

            % PROCESS AUTOCORRELATION
            ThreadCora = [];
            
            StackFilt = cat(3,StackFilt,zeros(Yt,Xt,Nfft-Zt));
            StackFilt = abs(ifft(abs(fft(StackFilt,[],3)).^2,[],3));
            Cora = squeeze(mean(mean(StackFilt,1),2));
            Cora = Cora(1:MaxSamp);
            Cora = Cora./(max(Cora));
            

            Correlations(Ind,:) = Cora;
            disp(horzcat('IIR Opt Processing: ',num2str(Counter),'/',num2str(Nf*No*Nc)));
            Counter = Counter + 1;
        end   
   end
end

% FIND OPTIMAL FILTER
Differences = zeros(1,Nf*No*Nc);
for i = 1:Nf*No*Nc
    
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
    
    Differences(i) = (10*PeakCoef + 1) * mean(Correlations(i,:));
   % disp(horzcat('Peak: ',num2str(10*PeakCoef),' Cf: ',num2str(TempCfs(i)),' Order: ',num2str(TempOrders(i))));

end


% F1 W1 W2 W3 ... F2 W1 W2 W3 ... F3 W1 W2 W3 ...
[~,I] = min(Differences);
OrderOpt = TempOrders(I);
FiltTypeOpt = TempFilts{I};
CfOpt = TempCfs(I);
IndexOpt = I;

if(CfOpt == Cfs(end))
    warning('High cut-off frequency autoselect');
end


end

