function [T1,T2,T3] = getUnmixCoefAlexa488XAlexa647Exp(pnc,fn)

    % so far only for dichroic ZT594RDC
    if strfind(pnc,'594')
        T1 = 0.03;
        T3 = 0.99;
        if strfind(fn,'Atto565')
            T2 = 0.44;
        elseif strfind(fn,'Flip565')
            T2 = 0.26; disp('Flip565: Experimental determination does not exist, use theoretical value!')
        elseif strfind(fn,'568')
            T2 = 0.57;
        elseif strfind(fn,'JF549')
            T2 = 0.37;
        end
    end
end
