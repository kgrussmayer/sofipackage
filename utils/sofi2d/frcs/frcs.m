% estimate resolution using sectorial Fourier ring correlation (sFRC)
% tomas.lukes@epfl.ch

function [resMid, resHigh, resLow, frcCurves] = frcs(in1,in2, settings)

fprintf('\n -- computing FRC curve--\n')

sc_in1 = max(in1(:));
sc_in2 = max(in2(:));
in1 = in1./sc_in1;
in2 = in2./sc_in2;

if settings.bcgsub > 0
    bcg1 = median(in1(:));
    bcg2 = median(in2(:));
    in1 = imadjust(in1,[settings.bcgsub*bcg1 1],[0 1]);
    in2 = imadjust(in2,[settings.bcgsub*bcg2 1],[0 1]);
end

sector = [0,1*pi/6;
    1*pi/6, 2*pi/6;
    2*pi/6, 3*pi/6;
    3*pi/6, 4*pi/6;
    4*pi/6, 5*pi/6;
    5*pi/6, 6*pi/6];

for isector = 1:size(sector,1)+1
    
    if isector==size(sector,1)+1
        % standard FRC over the whole circle
        %             [resolution, resolution_h, resolution_l,frc_curve] = frcsec(in1,in2,[]); % FRC using radsum in matlab
        [resolution, resolution_h, resolution_l,frc_curve] = frcfull(in1,in2); % original FRC
        disp('full circle');
        
    else
        % calculate FRC in the given sector
        [resolution, resolution_h, resolution_l,frc_curve] = frcsec(in1,in2,sector(isector,:));
        disp(['sector ',num2str(isector),'/',num2str(size(sector,1))]);
    end
    
    resMid(isector) = resolution*(settings.pixelsize);
    resHigh(isector) = resolution_h*(settings.pixelsize);
    resLow(isector)  = resolution_l*(settings.pixelsize);
    frcCurves{isector} = frc_curve;
end
% eof