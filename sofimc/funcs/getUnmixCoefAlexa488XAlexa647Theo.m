function [T1,T2,T3] = getUnmixCoefAlexa488XAlexa647Theo(pnc,fn)

    % so far only for dichroic ZT594RDC
    if strfind(pnc,'594')
        R1 = 0.02;
        T1 = 0.02;
        R3 = 57.45;
        T3 = 0.98;
        if strfind(fn,'Atto565')
            R2 = 0.55;
            T2 = 0.35;
        elseif strfind(fn,'Flip565')
            R2 = 0.34;
            T2 = 0.26;
        elseif strfind(fn,'568')
            R2 = 0.9;
            T2 = 0.47;
        elseif strfind(fn,'JF549')
            R2 = 0.19; 
            T2 = 0.16;
        end
    end
end
