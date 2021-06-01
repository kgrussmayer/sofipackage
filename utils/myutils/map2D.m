function im2D = map2D(im1D)


N = length(im1D);
x = linspace(-N/2,N/2,N);
[X,Y] = meshgrid(x);
r = sqrt(X.^2 + Y.^2);
r(r > max(x)) = max(x);
r = r+N/2;

rc  = ceil(r);
rf = floor(r);

p = rc-r;

tc = im1D(rc);
tf = im1D(rf);

% linear interpolation
im2D = p.*tf + (1-p).*tc;
end