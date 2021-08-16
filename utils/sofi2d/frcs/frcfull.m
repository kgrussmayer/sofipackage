function [resolution, resolution_h, resolution_l,frc_curve] = frcfull(in1,in2)
% calculate FRC in a given sector of the Frequency space
%
% inputs:
%
% in1 ...   first image (2D matrix)
% in2 ...   second image (2D matrix)
% sector ... angle for circular sum in rad, [start angle, end angle]
%            for example [0, pi/2]
%
% tomas.lukes@epfl.ch

in1 = dip_image(in1);
in2 = dip_image(in2);

sz = imsize(in1);

% Calculate FRC curve
% Compute mask in x-direction
nfac = 8;  % Image width / Width of edge region
x_im = xx(sz(1),sz(2))/sz(1);
mask = 0.5-0.5*cos(pi*nfac*x_im);
mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;

% Check that input images are square and mask
if sz(1) == sz(2)
    mask = mask*rot90(mask);
    
    % Mask input images
    in1 = mask*in1;
    in2 = mask*in2;
else
    warning('frc:nonsquare','Images are not square.');
    
    % Compute mask in y-direction
    y_im = yy(sz(1),sz(2))/sz(2);
    mask_y = 0.5-0.5*cos(pi*nfac*y_im);
    mask_y(abs(y_im)<((nfac-2)/(nfac*2))) = 1;
    mask = mask*mask_y;
    clear mask_y
    
    % Mask input images
    in1 = mask*in1;
    in2 = mask*in2;
    
    % Make images square through zero padding
    in1 = extend(in1,[max(sz) max(sz)]);
    in2 = extend(in2,[max(sz) max(sz)]);
end

% Fourier transform input images
in1 = ft(in1);
in2 = ft(in2);

% Compute fourier ring correlation curve
frc_num = real(radialsum(in1.*conj(in2)));              % Numerator
in1 = abs(in1).^2;
in2 = abs(in2).^2;
frc_denom = sqrt(abs(radialsum(in1).*radialsum(in2)));	% Denominator
frc_out = double(frc_num)./double(frc_denom);           % FRC
frc_out(isnan(frc_out)) = 0;                            % Remove NaNs

% Calculate the resolution
sz = max(imsize(in1));

% resolution = frctoresolution(frc_out,sz);
[resolution, resolution_h, resolution_l,frc_curve] = frctoresolution(frc_out,sz);
% eof