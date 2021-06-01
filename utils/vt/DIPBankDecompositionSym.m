function [ BankSignals, fig ] = DIPBankDecompositionSym( Stack, S, Nfilt, FiltOrder, Threads )
% SYMMETRIC BANK DECOMPOSITION
% Stack ...   Image Stack
% S .......   Input Parameters
% Nfilt ...   Number of prefered  filters
% FiltOrder . FIR Order used for decomposition
% Threads ... Number of working threads

[Y,X,Z] = size(Stack);
NumFiltLevels = floor(log2(Nfilt));
N = 2^(NumFiltLevels);
MaxSamp = S.MaxCorrSamp;

% CONSTRUCT LP & HP FIR FILTERS
al = 1;ah = 1;
bl = fir1(FiltOrder,0.5,'low');
bh = fir1(FiltOrder,0.5,'high');
fig = [];
BankSignals = zeros(N,Z/N,X*Y);

% LOOP ACROSS IMAGE DIMENSIONS
for i = 1:Y  % parfor(i = 1:Y,Threads)
    for j = 1:X
        Lindex = j + (i-1)*X;
        
        SignalsNew = [];
        % SUBBAND DECOMPOSITION
        PixSignals = squeeze(Stack(i,j,:))';
        for NFL = 1:NumFiltLevels


            % SAVE DATA
            if(NFL > 1)
                PixSignals = SignalsNew;
                SignalsNew = [];
            end

            % NUMBER OF SUBLEVELS
            NSub = 2^(NFL-1);
            %disp(horzcat('NumSub DownStep: ',num2str(NSub)));


            for ins = 1:NSub
                % LOAD SIGNAL
                Sig = PixSignals(ins,:);

                % PROCESS FILTER
                SigL = filter(bl,al,Sig);
                SigH = filter(bh,ah,Sig); 

                % DECIMATE
                SignalsNew(2*(ins-1) +1,:) = downsample(SigL,2);
                SignalsNew(2*ins,:) = downsample(SigH,2);

            end

        end  % DECOMPOSITION END
       
        % SAVE INTO ARRAY
        BankSignals(:,:,Lindex) = SignalsNew;
        
    end
end

[Nb,Nbs,NSigs] = size(BankSignals);
% Number of subbands
% Length of signal within bank
% Total number of signals

 SpecAxis = linspace(-1,1,Nbs);
 CoraAxis = 0:MaxSamp-1;
 
% SHOW SPECTRUM AND CORRELATIONS OF ALL THE BANKS
BankSpec = zeros(Nb,Nbs);
BankCor = zeros(Nb,Nbs * 2);

 parfor(i = 1:Nb,Threads)
     for j = 1:NSigs
         BankSig = squeeze(BankSignals(i,:,j));
         BankSpec(i,:) = BankSpec(i,:) + 20*log10(fftshift(abs(fft(BankSig))));
         BankCor(i,:) = BankCor(i,:) + abs(ifft( abs(fft([BankSig,zeros(1,Nbs)])).^2 ));
     end
 end
 % AVG COR AND SPEC
BankSpec = BankSpec./NSigs;
BankCor = BankCor(:,1:MaxSamp)./NSigs;
for i = 1:Nb
    BankCor(i,:) = BankCor(i,:)./max(BankCor(i,:));
end
 
 % PLOT FIGURE
 Cls = jet(Nb);
 fig = figure('name','Symmetric Bank Decomposition','position',[S.Width/2 - S.Wi/2,S.Height/2-S.He/2,S.Wi,S.He]);
 % SPECTRUM
 subplot(2,1,1);
 for  i = 1:Nb
      plot(SpecAxis,BankSpec(i,:),'color',Cls(i,:),'linewidth',S.Lw);
      hold on;grid on;
      Legends{i} = horzcat('Subband: ',num2str(i));
 end
xlabel('Normalized frequency x \pi rad/sample','FontSize',S.FontSize);
ylabel('20log_{10}|S(\omega)|','FontSize',S.FontSize);
legend(Legends);

% CORRELATIONS
subplot(2,1,2);
for i = 1:Nb
    plot(CoraAxis,BankCor(i,:),'color',Cls(i,:),'linewidth',S.Lw);
    hold on; grid on;
end
xlabel('Frame [-]','FontSize',S.FontSize);
ylabel('Correlation [-]','FontSize',S.FontSize);
legend(Legends);



end

