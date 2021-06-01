%Flatten raw cumulants - distance factor correction in all 3 spatial
%dimensions
%
%Inputs:
% sofi      Raw cumulants
% orders    Cumulant orders                     {all}
% stddev    Flatten standard deviations if set  {0:flatten mean}
%
%Output:
% sofi      Flat cumulants

function sofi=sofiAllFlatten3D(sofi,orders,stddev)
if nargin < 2 || isempty(orders)
   orders=find(~cellfun(@isempty,sofi));
end
orders=orders(orders > 1);
if isempty(orders)
   return;
end
if nargin < 3 || isempty(stddev) || ~stddev
   stddev=@mean;
else
   stddev=@std;
end
for order=orders(:).'
   image=(sofi{order});
   [X,Y,Z,N]=size(image);
   term=ones(order);
   for x=1:order
      for y=1:order
          for z = 1:order
            t=image(x:order:X,y:order:Y,z:order:Z,:);
            term(x,y,z)=feval(stddev,t(:));
          end
      end
   end
   term=mean(term(:))./term;
   term=repmat(term,ceil([X Y Z]/order));
   sofi{order}=image(:,:,:,:).*term(1:X,1:Y,1:Z,ones(1,N));
end
