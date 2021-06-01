function h = getAngularHistogram(mcSOFI2)

ch1 = mcSOFI2(:,:,1); % ch1(ch1 < 0.01) = [];
ch2 = mcSOFI2(:,:,2); % ch2(ch2 < 0.01) = [];
ch3 = mcSOFI2(:,:,3); % ch3(ch3 < 0.01) = [];

s1 = (ch2(:)./ch1(:));
s2 = (ch3(:)./ch2(:));
% s3 = (ch1(:)./ch3(:));

hx = linspace(-2,60,500);
h1 = hist(s1,hx);
h2 = hist(s2,hx);
% h3 = hist(s3,hx);

h1(1) = 0; h1(end) = 0;
h2(1) = 0; h2(end) = 0;
h = h1 + h2; %hist(s1+s2,hx);

figure(101);
plot(hx,h,'linewidth',2); hold on
plot(hx,h1,'linewidth',1.2);
plot(hx,h2,'linewidth',1.2);
% plot(hx,h3,'linewidth',1.2);
hold off
% xlim([-3.2 3.2])
legend('sum','ch2/ch1','ch3/ch2')
