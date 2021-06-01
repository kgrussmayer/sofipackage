function [ XX ] = DIPCalculation( Stack,StackFilt, S )
% CALCULATE SPECTRUM AND CORRELATION ON BOTH ORIGINAL AND CORRECTED DATA
% IN
% Stack ... Original Image Stack.
% StackFilt ... Corrected ImageStack via filter
% StackExp ... Corrected ImageStack via inv. fit
% S ... All options and settings
%
% OUT
% Spec = Average Complex Spectrum
% Phase = Average Phase of the Spectrum in Degrees
% Corela = Average Correlation in Time
% SpecFilt = Average Spectrum of the corrected Image Data
% PhaseFilt = Average Phase of the corrected Image Data 
% CorelaFilt = Average Correlation of the corrected Image Data
% SpecExp = Average Spectrum of the corrected Image Data
% PhaseExp = Average Phase of the corrected Image Data 
% CorelaExp = Average Correlation of the corrected Image Data


% EXTRACT PARAMS
[Y,X,Z] = size(Stack);
% Wsize = S.WSize;
MaxCorrSamp = S.MaxCorrSamp;


% CALC OPTIONS
SC = S.SC;
SCSpec = S.SCSpec;
SCPhase = S.SCPhase;
% SCCorr = S.SCCorr;


XX.Spec = [];
XX.SpecFilt = [];
XX.SpecExp = [];
XX.Phase = [];
XX.PhaseFilt = [];
XX.PhaseExp = [];
XX.Corela = [];
XX.CorelaFilt = [];
XX.CorelaExp = [];

% SPECTRUM
if(SCSpec && SC)
    XX.SpecFilt = squeeze(mean(mean(fft(StackFilt,[],3),2),1));
%     XX.SpecExp  = squeeze(mean(mean(fft(StackExp,[],3),2),1));
    XX.Spec     = squeeze(mean(mean(fft(Stack,[],3),2),1));
end

% PHASE SPECTRUM
if(SCPhase && SC)
    XX.Phase =      squeeze(mean(mean( rad2deg(angle( fft(Stack,[],3)))     ,2),1));
    XX.PhaseFilt =  squeeze(mean(mean( rad2deg(angle( fft(StackFilt,[],3))) ,2),1));
%     XX.PhaseExp =   squeeze(mean(mean( rad2deg(angle( fft(StackExp,[],3)))  ,2),1));
end


    CorelaFilt = abs(ifft(abs(fft(cat(3,StackFilt,zeros(Y,X,Z)),[],3)).^2,[],3));
    CorelaFilt = squeeze(mean(mean(CorelaFilt,1),2));
    XX.CorelaFilt = CorelaFilt(1:MaxCorrSamp);
    
    Corela = abs(ifft(abs(fft(cat(3,Stack,zeros(Y,X,Z)),[],3)).^2,[],3));
    Corela = squeeze(mean(mean(Corela,1),2));
    XX.Corela = Corela(1:MaxCorrSamp);
%     
%     CorelaExp = abs(ifft(abs(fft(cat(3,StackExp,zeros(Y,X,Z)),[],3)).^2,[],3));
%     CorelaExp = squeeze(mean(mean(CorelaExp,1),2));
%     XX.CorelaExp = CorelaExp(1:MaxCorrSamp);


end

