function [ BankSignals, fig ] = DIPBankDecompositionASym( Stack, S, Nfilt, FiltOrder, Threads )
% ASYMMETRIC BANK DECOMPOSITION
% Stack ...   Image Stack
% S .......   Input Parameters
% Nfilt ...   Number of prefered  filters
% FiltOrder . FIR Order used for decomposition
% Threads ... Number of working threads

[Y,X,Z] = size(Stack);
MaxSamp = S.MaxCorrSamp;

% CONSTRUCT LP & HP FIR FILTERS
al = 1;ah = 1;
bl = fir1(FiltOrder,0.5,'low');
bh = fir1(FiltOrder,0.5,'high');
BankSignals = cell(Y,X);
fig = [];

% LOOP ACROSS IMAGE DIMENSIONS
parfor (i = 1:Y,Threads)  % parfor(i = 1:Y,Threads)
    for j = 1:X
        
        
        % SUBBAND DECOMPOSITION
        SignalsNew = []; PixSignals = [];
        PixSignals{1} = squeeze(Stack(i,j,:))';
        for NFL = 1:Nfilt


            % SAVE DATA
            if(NFL > 1)
                PixSignals = SignalsNew;
                SignalsNew = [];
            end

            Ns = length(PixSignals);
            
            % LOAD SIGNAL
            Sig = PixSignals{1};
            
            % PROCESS FILTER
            SigL = filter(bl,al,Sig);
            SigH = filter(bh,ah,Sig);
            
            % DECIMATE
            SignalsNew{1} = downsample(SigL,2);
            SignalsNew{2} = downsample(SigH,2);
            
            % COPY THE REST
            for ii = 1:Ns -1
                SignalsNew{ii+2} = PixSignals{ii+1};
            end


        end  % DECOMPOSITION END
       
        % SAVE INTO ARRAY
        BankSignals{i,j} = SignalsNew;
        
        
    end
end

NSigs = length(BankSignals);
Nb = length(BankSignals{1,1});

% GET LENGTHS OF THE SPECTRA
T = BankSignals{1,1};
parfor i = 1:Nb
    Lengths(i) = length(T{i});
end



 
% SHOW SPECTRUM AND CORRELATIONS OF ALL THE BANKS
BankSpec = cell(Nb);
BankCor = cell(Nb);

% Preinit spec and corr Cells
parfor i = 1:Nb
    BankSpec{i} = zeros(1,Lengths(i));
    BankCor{i} = zeros(1,2*Lengths(i));
end

for i = 1:Y
    for j = 1:X
        BankSigs = BankSignals{i,j};
        
        for sub = 1:Nb
          Sig = BankSigs{sub};
          %disp(horzcat('i: ',num2str(i),' j: ',num2str(j),' Sig: ',num2str(length(Sig)),' BankSpec: ',num2str(length(BankSpec{sub}))));
          BankSpec{sub} = BankSpec{sub} +  20*log10(fftshift(abs(fft(Sig))));
          BankCor{sub} = BankCor{sub} + abs(ifft( abs(fft([Sig,zeros(1,Lengths(sub))])).^2 ));
        end
    end
    
end

 % AVG COR AND SPEC
parfor (i = 1:Nb,Threads)
    BankSpec{i} = BankSpec{i}./NSigs;
    Cora = BankCor{i};
    Cora = Cora(1:min(length(Cora),MaxSamp));
    Cora = Cora./NSigs;
    Cora = Cora./max(Cora)
    BankCor{i} = Cora;
end
 
 % PLOT FIGURE
 Cls = jet(Nb);
 fig = figure('name','Symmetric Bank Decomposition','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
 
 % SPECTRUM
 subplot(2,1,1);
 for  i = 1:Nb
      SpecAxis = linspace(-1,1,Lengths(i));
      plot(SpecAxis,BankSpec{i},'color',Cls(i,:),'linewidth',S.Lw);
      hold on;grid on;
      Legends{i} = horzcat('Subband: ',num2str(i));
 end
xlabel('Normalized frequency x \pi rad/sample','FontSize',S.FontSize);
ylabel('20log_{10}|S(\omega)|','FontSize',S.FontSize);
legend(Legends);

% CORRELATIONS
subplot(2,1,2);
for i = 1:Nb
    CoraAxis = 0:length(BankCor{i}) - 1;
    plot(CoraAxis,BankCor{i},'color',Cls(i,:),'linewidth',S.Lw);
    hold on; grid on;
end
xlabel('Frame [-]','FontSize',S.FontSize);
ylabel('Correlation [-]','FontSize',S.FontSize);
legend(Legends);



end

