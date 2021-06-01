function [ b,a,ht,hf,Nt,Fig ] = DIPFilterConstruction( FType,Type,Cf1,Cf2,Ord,PassRipple,StopRipple,ShowFilt,FirWind )
% FILTER CONSTRUCTION FUNCTION
% IN:
% FType = Filter Frequency Type (High/Low/Pass/Stop)
% Type = Filter Type (FIR/Butter/Cheby/InvCheby/Cauer)
% Cf1 = Normalized CutOffFrequency [0->1] (All types)
% Cf2 = Normalized CutOffFrequency [0->1] (Only Pass/Stop)
% Ord = Filter Order [-]
% PassRipple [dB] (Chebyshev /Cauer)
% StopRipple [dB] (Inv Chebyshev / Cauer)
% ShowFilt = To Show Filter Response
% FirWind - Optional parameter for fir filter design
%
% OUT:
% b - Filter Coeficients
% a - Filter Coeficients
% ht - Impulse Response
% hf - Frequency Response
% Nt - Number of samples in the Impulse Response 

if(mod(Ord,2)== 0)
    NWin = Ord+1;
else
    NWin = Ord+2;
end

% Fill Unused arguments
switch nargin
    case 5
        PassRipple = 0;
        StopRipple = 0;
        ShowFilt = 0;
        FirWind = nuttallwin(NWin)';
    case 6
        StopRipple = 0;
        ShowFilt = 0;
        FirWind = nuttallwin(NWin)';
    case 7
        ShowFilt = 0;
        FirWind = nuttallwin(NWin)';
    case 8
        FirWind = nuttallwin(NWin)';
end


Cl1 = [127,0,204]/255;     % Violet
Cl2 = [0,128,255]/255;     % Blue

% ECC CHECK
if(Cf1 > 1 || Cf1 < 0 )
    error('Please provide Normalized Frequency for Filter Design!');
end

%disp(horzcat('Designing ',FType,' Filter'));



if(strcmp(FType,'high') || strcmp(FType,'low'))
    % HIGH PASS OR LOWPASS
    Tit =[]; 
    switch Type
        case 'fir'
            a = 1;b = fir1(Ord,Cf1,FType,FirWind);
            Tit = 'FIR';
        case 'butter'
            [b,a] = butter(Ord,Cf1,FType);
            Tit = 'ButterWorth';
        case 'cheby'
            [b,a] = cheby1(Ord,PassRipple,Cf1,FType);
            Tit = 'Chebyshev';
        case 'invcheby'
            [b,a] = cheby2(Ord,StopRipple,Cf1,FType);
            Tit = 'Inverse Chebyshev';
        case 'cauer'
            [b,a] = ellip(Ord,PassRipple,StopRipple,Cf1,FType);
            Tit = 'Cauer';
        otherwise
            error(horzcat('Wrong Filter! Received: ',Type));
    end
elseif (strcmp(FType,'pass') || strcmp(FType,'stop'))
    if(strcmp(FType,'pass'))
        FType = 'bandpass';
    end
    % BAND PASS OR BAND STOP
    Tit = [];
    switch Type
        case 'fir'
            a = 1;b = fir1(Ord,[Cf1,Cf2],FType,FirWind);
            Tit = 'FIR';
        case 'butter'
            [b,a] = butter(Ord,[Cf1,Cf2],FType);
            Tit = 'ButterWorth';
        case 'cheby'
            [b,a] = cheby1(Ord,PassRipple,[Cf1,Cf2],FType);
            Tit = 'Chebyshev';
        case 'icheby'
            [b,a] = cheby2(Ord,StopRipple,[Cf1,Cf2],FType);
            Tit = 'Inverse Chebyshev';
        case 'cauer'
            [b,a] = ellip(Ord,PassRipple,StopRipple,[Cf1,Cf2],FType);
            Tit = 'Cauer';
        otherwise
            error(horzcat('Wrong Filter! Received: ',Type));
    end
else
    error('Unknown Filter Type! Expecting High/Low/Pass/Stop')
end

% CALCULATE FREQ. RESPONSE
[hf,w] = freqz(b,a);
[ht,t] = impz(b,a);
Nt = length(ht);
if(Nt > 3000)
    warning('Long filter impulse response - May Have a negative effect on the image');
end

% SHOW FILTER
if(isequal(ShowFilt,1))
    Fig = figure('name','Frequency response of the selected filter');
    subplot(2,1,1);
    plot(w/pi,20*log10(abs(hf)),'color',Cl1,'linewidth',1);
    xlabel('Normalized Frequency (x \pi rad/sample)');
    ylabel('20*log10(|H(\omega)|)');
    title(horzcat(Tit,' ',num2str(Ord),' Order Frequency Response'));
    grid on;
    subplot(2,1,2);
    stem(t,ht,'color',Cl2,'marker','.','markersize',17);
    xlabel('Delay [Samples]');
    ylabel('|h(\tau)|');
    title(horzcat(Tit,' ',num2str(Ord),' Order Impulse Response'));
else
    Fig = [];
end


end

