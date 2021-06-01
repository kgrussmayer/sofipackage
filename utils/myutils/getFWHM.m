function sx = getFWHM(dat,display)

if nargin < 2; display = 0; end
% return the exact FWHM of dat using local linear interpolation
% dat has to be a vector of any arbitrary focus


try
%     dat = linmap(dat,min(dat),max(dat),0,1);
	% substract background
dat = dat - min(dat);

M = max(dat);
y = M/2;

map = abs(diff(dat > y));
temp = find(map == 1);

xa1 = temp(1);xb1 = temp(2);
xa2 = xa1+1; xb2 = xb1+1;

ya1 = dat(xa1); ya2=dat(xa2);
yb1=dat(xb1); yb2=dat(xb2);

ma = ya2-ya1;mb = yb2-yb1;
ca = ya1-ma*xa1;cb = yb1-mb*xb1;

sx = abs((y-ca)/ma - (y-cb)/mb);

% check the analysis by plot
if display
    figure(display)
    plot(dat)
    hold on
    plot(dat,'*')
    plot([1 length(dat)],[y y])
    plot(xa1,ya1,'g*');plot(xa2,ya2,'g*','linewidth',2)
    plot(xb1,yb1,'g*');plot(xb2,yb2,'g*','linewidth',2)
    plot((y-ca)/ma,y,'r*');plot((y-cb)/mb,y,'r*')

%     xlim([length(dat)/2-40 length(dat)/2+40]);ylim([0.4 0.6])

    hold off
end

catch
    sx = length(dat);
end