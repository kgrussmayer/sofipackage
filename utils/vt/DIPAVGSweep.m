function [ OrdOptimum,Correlations] = DIPAVGSweep( StackO, S ,Orders, Px, Threads )

% FIND MAXIMUM
MaxSamp = S.MaxCorrSamp;
Widefield = sum(StackO,3);
[M,I] = max(Widefield);
[MM,II] = max(M);
Yi = I(II);
Xi = II;
N = length(Orders);

% FIX EDGES PROBLEMS
[Y,X,Z] = size(StackO);
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

for k  = 1:N
    disp(horzcat('Moving Average Sweep: ',num2str(k),'/',num2str(N)))
    Ord = Orders(k);
    af = 1; bf = 1/Ord*ones(1,Ord);
    Ht = impz(bf,af);
    Nh = length(Ht);
    Nfft = 2^(ceil(log2(Nh + Z - 1)));
    Hzp = fft([Ht',zeros(1,Nfft - Nh)]);
    

    % FILTER STACK
     [Yt,Xt,Zt] = size(StackC);
%     StackC = cat(3,StackC,zeros(Yt,Xt,Nfft-Zt));
%     StackFilt = fft(StackC,[],3);
%     for i = 1:Nfft
%         StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
%     end
%     StackFilt = real(ifft(StackFilt,[],3));
%     StackC = StackC(:,:,1:Zt);
    StackFilt = StackC;
    parfor (i = 1:Yt,Threads)
        for j = 1:Xt
            OSignal = squeeze(StackC(i,j,:))';
            MeanSig = filter(bf,af,OSignal);
            MeanSig = circshift(MeanSig,[0,-Ord/2]);
            StackFilt(i,j,:) = OSignal - MeanSig;
        end
    end
    
    % CUT AFTER FILTRATION
    StackFilt = StackFilt(:,:,S.ImStart-Ord/2:S.ImStart + S.NumImages-Ord/2 - 1);
    [Y,X,Z] = size(StackFilt);

    % COMPUTE CORRELATION
    CorelaFilt = cat(3,StackFilt,zeros(Y,X,Z));
    CorelaFilt = abs(ifft(abs(fft(CorelaFilt,[],3)).^2,[],3));
    CorelaFilt = squeeze(mean(mean(CorelaFilt,1),2));
    CorelaFilt = CorelaFilt(1:S.MaxCorrSamp);
    CorelaFilt = CorelaFilt./max(CorelaFilt);
    Correlations(k,:) = CorelaFilt;
    
end

%% FIND OPTIMAL FILTER BASED ON ORDER
Differences = zeros(1,N);
for i = 1:N    
    
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

[~,I] = min(Differences);
OrdOptimum = Orders(I);


end

