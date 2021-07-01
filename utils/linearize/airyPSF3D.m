function psf = airyPSF3D(fwhm)
% return an ideal 3D psf model where fwhm is a vector
% fwhm = [fwhmx, fwhmy, fwhmz];

kcx = 1.815/fwhm(1);
kcy = 1.815/fwhm(2);
kcz = 1.6/fwhm(3);

Nx = ceil(3*fwhm(1)); Nx = Nx + not(mod(Nx,2)); Nx = max(Nx,15);
Ny = ceil(3*fwhm(2)); Ny = Ny + not(mod(Ny,2)); Ny = max(Ny,15);
Nz = ceil(3*fwhm(3)); Nz = Nz + not(mod(Nz,2)); Nz = max(Nz,11);
upsc = 4;
x = linspace(-1,1,upsc*Nx); y = linspace(-1,1,upsc*Ny); z = linspace(-1,1,upsc*Nz);
[X,Y,Z] = ndgrid(x,y,z); R = sqrt(X.^2 + Y.^2);

% CTF approximation using a sine function 
kzmask = kcz.*sin((pi/kcx).*R);
ctf = (abs(Z)-0.01 < kzmask).*(R < kcx);

% ctf balancing such that its projection is an ideal 2D psf
R = R(:,:,1);
weight = ((1-R./kcx)./sum(ctf,3));
ctf = ctf.*repmat(weight,[1 1 size(ctf,3)]);
ctf(isnan(ctf)) = 0;

% psf is the inverse Fourier transform of the CTF
psf = linmap(abs(ifftshift(ifftn(ctf))),0,1);

% center the psf
psf(1,:,:) = []; psf(:,1,:) = [];  psf(:,:,1) = []; 
psf(end+1,:,:) = psf(1,:,:);
psf(:,end+1,:) = psf(:,1,:);
psf(:,:,end+1) = psf(:,:,1);
Nx = (Nx-1)/2; Ny = (Ny-1)/2; Nz = (Nz-1)/2;
psf = psf(end/2-Nx:end/2+Nx,end/2-Ny:end/2+Ny,end/2-Nz:end/2+Nz);
% eof