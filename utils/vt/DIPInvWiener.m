function [ StackOut ] = DIPInvWiener(IM,PSF )
% INVERSE WIENER FILTRATION

% Norm
IMX = IM./max(IM(:));

[Hy,Hx] = size(PSF);
[YY,XX] = size(IMX);
Nx = XX + 2*Hx;
Ny = YY + 2*Hy;

% CREATE OTF
OTF = psf2otf(PSF,[Ny,Nx]);

 % ESTIMATE NOISE
Sig = mean(IMX(:));
Noise = std(IMX(:));
SNR = 20*log10(Sig/Noise);
% if(SNR<20)
%     SNR = 20;
% end


% ADD BORDERS
IMX = padarray(IMX,[Hy,Hx],'symmetric'); 
IMX = fft2(IMX);

% WIENER
H = conj(OTF) ./ (abs(OTF).^2 + 1/10^(SNR/10));
IMX =  real((ifft2(IMX.*H)));
IMX = IMX + min(IMX(:));

% REMOVE BORDERS
IMX = IMX(Hy:YY+Hy-1,Hx:XX+Hx-1);
   
% SAVE
StackOut = abs(IMX);

end

