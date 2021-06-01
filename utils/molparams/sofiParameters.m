function [ratio,density,bright]=sofiParameters(sofi,tirf)
%Estimate emitter parameters from flat cumulants.
%
%Inputs:
% sofi  ...  cell  Flat cumulants of orders 2 to 4
% tirf  ...  bool  TIRF illumination if true % TODO: needs better explanation!
%
%Outputs:
% ratio     On-time ratio
% density   Emitter density
% bright    Emitter brightness

if nargin > 1 && tirf
   tirf=2;
else
   tirf=1.5;
end
if numel(sofi) < 4 || any(cellfun(@isempty,sofi(2:4)))
   error('sofi:orders','Require flat cumulant images of 2nd, 3rd and 4th order.');
end
cum2=(sofi{2}(:,:,:));
cum3=(sofi{3}(:,:,:));
cum4=(sofi{4}(:,:,:));
%
% Resample cumulants to highest order.
%
[x,y]=xygrid(size(cum4));
cum2=interpolate(cum2,x,y);
cum3=interpolate(cum3,x,y);
%
% Estimate the parameters.
%
k1=1.5^tirf*cum3./cum2;
k2=2^tirf*cum4./cum2;

t0=1./sum(cum2,3);

t1=sqrt(max(0,3*k1.^2 - 2*k2));

bright=sum(cum2.*t1,3).*t0;

bright(bright<0) = 0;

t2=max(0,min(0.5-0.5*k1./t1,1));

ratio=sum(cum2.*t2,3).*t0;
ratio(ratio<0)=0;
ratio(ratio>1)=1;

t2=max(0,cum2./(t1.^2.*t2.*(1-t2)));
t2(isinf(t2)) = 0;

density=sum(cum2.*t2,3).*t0;

% ratio2 = 3*k1.^2 - k1.*sqrt(3*k1.^2 + 2*k2) - 2*k2;
% ratio2 = ratio2./(2*(3*k1.^2-2*k2));

%%
% for n = 1:size(sofi,3)
%     k1=sofi{3}(1:3:end,1:3:end,n)./sofi{2}(1:2:end,1:2:end,n);
% %     k1=interpolate(sofi{3}(:,:,n),x,y)./interpolate(sofi{2}(:,:,n),x,y);
%     
%     k2=sofi{4}(1:4:end,1:4:end,n)./sofi{2}(1:2:end,1:2:end,n);
% % k2=sofi{4}(:,:,n)./interpolate(sofi{2}(:,:,n),x,y);
% 
%     A=(3*k1.^2-2*k2);
%     A(A<=0)=0;
%     A=sqrt(A);
% 
%     rate = 0.5-0.5*k1./A;
%     rate(rate<0)=0;
%     rate(rate>1)=1;
% 
%     Nmol=sofi{2}(1:2:end,1:2:end,n)./A.^2./(rate.*(1-rate));
% %     Nmol=interpolate(sofi{2}(:,:,n),x,y)./A.^2./(rate.*(1-rate));
%     Nmol(isinf(Nmol)) = 0;
% %     imagesc(Nmol);
%     Nmols(:,:,n) = Nmol;
% end
% 
% 
% imagesc(sum(abs(Nmols),3))
% temp = sum(abs(Nmols),3);
% sum(temp(:))

