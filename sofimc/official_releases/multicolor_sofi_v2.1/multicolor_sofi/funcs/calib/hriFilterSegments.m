function out=hriFilterSegments(Id,aupl,alol,sys,out)
% Description
% ===========
% hriFilterSegments chooses the segments that show the right size and
% eliminates accumulations and weak molecules. Further it computes
% the weighted center of gravity of each segment, which will deliver an
% estimation on the position of the molecule.
%
% Input
% =====
% Id:   image to be segmented
% aupl: Upper limit of the segments area
% alol: Lower limit of the segments area
% sys:  The system variables
% out:  The output variables
%
% Output
% ======
% out.xoAbsolute,out.yoAbsolute:    The absolute positions of the weighted
%                                   center of gravity of the segments in 
%                                   pixels
% out.xoEstimate,out.yoEstimate:    The estimated deviations from the
%                                   absolute positions [1/10 pixels]
% out.d: The segments (images of single molecules)
%
% Author
% ======
%      Stefan Geissbuehler
%      Swiss Federal Institute of Technology, CH-1015 Lausanne
%      Laboratoire d'optique biomedicale, LOB
%      Biomedical imaging group, BIG

%filter out
a=find(out.segA<aupl);
Af=out.segA(a);
xf=out.segx(a);
yf=out.segy(a);
a=find(Af>alol);
xf=xf(a);
yf=yf(a);

a=out.ru-1+out.nh;

b=(xf>a-1);
xf=xf(b);
yf=yf(b);
b=(xf<size(Id,1)-a+1);
xf=xf(b);
yf=yf(b);
b=(yf>a-1);
xf=xf(b);
yf=yf(b);
b=(yf<-(a-1)+size(Id,2));
xf=xf(b);
yf=yf(b);

out.xoAbsolute=round(xf);
out.yoAbsolute=round(yf);

out.xoEstimate=10*(xf-out.xoAbsolute);
out.yoEstimate=10*(yf-out.yoAbsolute);

out.x=(-a:a);
out.y=(-a:a);

out.d=sys.bg*ones(numel(out.x),numel(out.y),numel(out.xoEstimate));

for m=1:numel(out.xoEstimate)
    for n=1:numel(out.x)
        for o=1:numel(out.y)
            if(out.xoAbsolute(m)+out.x(n)>0 && out.xoAbsolute(m)+out.x(n)<size(Id,1) && out.yoAbsolute(m)+out.y(o)>0 && out.yoAbsolute(m)+out.y(o)<size(Id,2))
                out.d(n,o,m)=Id(out.xoAbsolute(m)+out.x(n),out.yoAbsolute(m)+out.y(o));
            end
        end
    end
end