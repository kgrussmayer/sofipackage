function [ FinSofi,Correlations ] = DIPSOFIBankFilter( Stack, S )
% PROCESS BANK FILTERING USING FIR FILTERS OF SPECIFIED ORDER, CFs and WINDOW FUCTION

ImStart = S.ImStart;
NumImg = S.NumImages;
Ord = S.FiltOrder;
CF1s = S.CF1;
CF2s = S.CF2;
MaxCorr = S.MaxCorrSamp;


% PRE - INIT
NumFilt = length(CF1s);
[Y,X,Z] = size(Stack);Ns = Z;


try
    % RELEASE WAITBAR ...
    WBar = waitbar(0,'1','Name','SOFI Computation in progress, please wait ...');

    for k = 1:NumFilt
        % CALCULATE FILTER COEF
        if(k == 1)
            % SIMPLE LOWPASS
            [b,a,Ht,~,Nh,~] = DIPFilterConstruction('low','fir',CF2s(k),[],Ord,[],[],[]);
        
        elseif(k == NumFilt)
            % SIMPLE HIGHPASS
            [b,a,Ht,~,Nh,~] = DIPFilterConstruction('high','fir',CF1s(k),[],Ord,[],[],[]);
        else
            % ANYTHING BETWEEN IS BANDPASS
            [b,a,Ht,~,Nh,~] = DIPFilterConstruction('pass','fir',CF1s(k),CF2s(k),Ord,[],[],[]);
        end
        
        % FILTER DATA
        Nfft = Ns + 2*Nh + 1;
        Hzp = fft([Ht',zeros(1,Nfft - Nh)]);
        StackFilt = fft(cat(3,Stack,zeros(Y,X,Nfft-Z)),[],3);
        for i = 1:Nfft
            StackFilt(:,:,i) = StackFilt(:,:,i)*Hzp(i);
        end
        StackFilt = real(ifft(StackFilt,[],3));
        StackFilt = StackFilt(:,:,ImStart:ImStart + NumImg - 1);
        
        % PROCESS AUTOCORRELATION
        StackCora = abs(fft(cat(3,StackFilt,zeros(Y,X,NumImg)),[],3));
        StackCora = squeeze(mean(mean(abs(ifft(StackCora.^2,[],3)),1),2));
        StackCora = StackCora(1:MaxCorr);
        Correlations(k,:) = StackCora./max(StackCora);

        % CALCULATE SOFI AND POSTPROCESS
        [SOFI,fwhm] = DIPSOFICumulants(StackFilt,S);
        [SOFI] = sofiLinearize(SOFI,fwhm,[],S.Orders,S.NumIter,S.Deconv);
        
        % ALL SOFI ORDERS
        for s = 1:length(S.Orders)
            Sofies{k,s} = SOFI{s};
        end
        
        % SHOW STATUS (INTEGRATION STATUS)
        waitbar(k / NumFilt,WBar,sprintf('%d//%d',k,NumFilt));
    end
    
    % REMOVE WAITBAR
    delete(WBar);
    
catch ME
    
    % CLOSE WAITBAR AND RETHROW
    if(exist('WBar','var'))
        delete(WBar);
    end
    rethrow(ME);
end


% JOIN ALL RESULTS
SSum = sum(Correlations,2);
[Ys,Xs,~] = size(Sofies);


% CALCULATE WEIGHTS ( The Higher Normalized AutoCorr, the lower weight)
Coefs = 1./SSum;
Coefs = Coefs/(sum(Coefs));
disp(horzcat('Filter bank Weights: ',num2str(Coefs')));

% MERGE RESULTS FROM BANKS AND DIFERENT ORDERS
for i = 1:NumFilt
    for s = 1:length(S.Orders)
        if(i == 1)
            FinSofi{s} = Coefs(i)*Sofies{i,s};
        else
            FinSofi{s} = FinSofi{s} + Coefs(i)*Sofies{i,s};
        end
    end
end



end

