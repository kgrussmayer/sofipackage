function [ Correlations, FilterTypes, FilterOrder ] = DIPAllFilter(StackO, S )
% FUNCTION FOR CALCULATING CORRELATIONS FOR ALL TYPES OF FILTERS
% - BUTTERWORTH
% - CHEBYSHEV
% - INVERSE CHEBYSHEV
% - CAUER
% - FIR

OrdFir = S.FiltOrderFir;
OrdIir = S.FiltOrderIir;
Cf = S.FiltCutOff1;
PassRip = S.FiltPassRip;
StopRip = S.FiltStopRip;
MaxSamp = S.MaxCorrSamp;
ImStart = S.ImStart;
NumImg = S.NumImages;
[Y,X,Z] = size(StackO);

Correlations = [];
FilterTypes = {'Fir','Butter','Cheby','Invcheby','Cauer'};
FilterOrder = [OrdFir,OrdIir*ones(1,4)];
Ns = length(StackO(1,1,:)); 

for k = 1:5
    disp(horzcat('Processing: ',FilterTypes{k},' - ',num2str(k),' / ',num2str(5)));
    
    %% FILTERING 
    [~,~,Ht,~,Nh,~] = DIPFilterConstruction('high',lower(FilterTypes{k}),Cf,0,FilterOrder(k),PassRip,StopRip,0);
    Nfft = 2^(ceil(log2(Nh + Z - 1)));
    Hzp = fft([Ht',zeros(1,Nfft-Nh)]);
    
    StackFilt = fft(cat(3,StackO,zeros(Y,X,Nfft - Z)),[],3);
    for  i =1:Nfft
        StackFilt(:,:,i) = StackFilt(:,:,i) * Hzp(i);
    end
    StackFilt = real(ifft(StackFilt,[],3));
    StackFilt = StackFilt(:,:,S.ImStart:S.ImStart + S.NumImages - 1);

    % PROCESS AUTOCORRELATION
    StackFilt = abs(fft(cat(3,StackFilt,zeros(Y,X,NumImg)),[],3));
    StackFilt = squeeze(mean(mean(abs(ifft(StackFilt.^2,[],3)),1),2));
    Correlations(k,:) = StackFilt./max(StackFilt);
    
end

% CUT CORRELATIONS
Correlations = Correlations(:,1:MaxSamp);

end

