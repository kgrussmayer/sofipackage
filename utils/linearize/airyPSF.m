function psf = airyPSF(fwhm)

kc = 1.815/fwhm;
N = ceil(8*fwhm); N = N + not(mod(N,2));
N = max(N,15);
x = linspace(-1,1,N);
ctf = (1-abs(x)./kc).*(abs(x)<kc);
psf = linmap(abs(ifftshift(ifft(ctf))),0,1);
psf(1) = []; psf(end+1) = psf(1);
% eof