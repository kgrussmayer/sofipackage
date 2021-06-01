function out = filterImCirc(in,r1,r2)

Nx = min(size(in,1),size(in,2));

[x,y] = meshgrid(linspace(-1,1,Nx));
r = sqrt(x.^2 + y.^2);
map = (r - r1)./(r2-r1).*(r < r2);
map(isnan(map)) = 0;
d = linmap(map,1,0,-pi/2,pi/2);
mask = (r < r1) + (sin(d)+1)/2.*(r < r2 & r >= r1);

out = real(ifft2(ifftshift(fftshift(fft2(in)).*mask)));