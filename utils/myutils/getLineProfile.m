function [trc,subIM] = getLineProfile(input,pa,pb,bin,display)


if nargin < 5; display = 0; end
if nargin < 4; bin = 4; end

off = round(sqrt(2)*length(input)/4);
input = padarray(input,[off off],0,'both');


pa = pa + off; pb = pb + off;
if (pa(1) < pb(1))
    p1 = pa; p2 = pb;
else
    p1 = pb; p2 = pa;
end
th = rad2deg(atan2(p2(2)-p1(2),p2(1)-p1(1)));

% rotate image
im = imrotate(input,(th),'bilinear','crop');

[m,n] = size(im);


% rotate points
M = [cosd(-th) -sind(-th) ; sind(-th) cosd(-th)];
c = [m/2 ; n/2];

p1 = M*(p1-c) + c;
p2 = M*(p2-c) + c;

x_start = min(round(p2(1)),round(p1(1))); x_stop = max(round(p2(1)),round(p1(1)));
y = round((p2(2)+p1(2))/2);

% clamp values within FOV
% x_start = clamp(x_start,1,n); x_stop = clamp(x_stop,1,n);
% y = clamp(y,1+bin,m-bin);

% sum over a the rows
subIM = im(y-bin:y+bin,x_start:x_stop);
trc = sum(subIM,1);

if display
    figure(display)
subplot(221);imagesc(input); hold on
    plot(pa(1),pa(2),'xr');
    plot(pb(1),pb(2),'xr');
    title('Input data'); hold off
subplot(222);imagesc(im);hold on
    plot(p2(1),p2(2),'xr');
    plot(p1(1),p1(2),'xr');
    rectangle('position',[x_start y-bin  ...
    x_stop-x_start 2*bin+1],'EdgeColor','w')
    title('Rotated image'); hold off
subplot(223);plot(trc)
end