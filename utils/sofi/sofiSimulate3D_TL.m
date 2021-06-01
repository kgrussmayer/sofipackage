%stack=sofiSimulate3D(frames,planes,centers,Ion,Ton,Toff,fwhm)
%-------------------------------------------------------------
%
%Simulate the acquisition of an image sequence of blinking emitters.
%
%Inputs:
% frames    Number of images
% planes    Number of focal planes
% centers   Emitter coordinates [x y z] [pixel]
% Ion       Signal per frame (on state) [photon]
% Ton       Average duration of the on state [frame]
% fwhm      PSF full width at half-maximum diameters [xy z] [pixel]
%
%Output:
% stack     Image sequence [101 x 101 x planes x frames]
%
function stack=sofiSimulate3D_TL(frames,centers,Ion,Ton,Toff,Tbl,sx,sy,sz,psf,b,gnoise)

indx = centers(:,1);
indy = centers(:,2);
indz = centers(:,3);

emitters=size(centers,1);

Ion=repmat(Ion(:),emitters/numel(Ion),1);

Ton=repmat(Ton(:),emitters/numel(Ton),1);

Toff=repmat(Toff(:),emitters/numel(Toff),1);

fig=statusbar('Simulation 3D: Photons...');

photons=zeros(emitters,frames);

stack=zeros(sx,sy,sz,frames);

for emitter=1:emitters
   photons(emitter,:)=brightness(Ion(emitter),Ton(emitter),Toff(emitter),Tbl,frames);
   stack(indy(emitter),indx(emitter),indz(emitter),:)=brightness(Ion(emitter),Ton(emitter),Toff(emitter),Tbl,frames);
   fig=statusbar(emitter/emitters,fig);
   
%    if isempty(fig);
%       return;
%    end
end

% x=cell(planes,1);
% y=cell(planes,1);
% n=-50:50;
% for plane=1:planes
%    w=fwhm(1).*sqrt(1 + 4/fwhm(2).^2.*(plane-(1+planes)/2 - centers(:,3)).^2);
%    x{plane}=gaussian(w,centers(:,1),n);
%    y{plane}=gaussian(w,centers(:,2),n).';
% end
fig=statusbar('Simulation 3D: Images...',fig);
% n=numel(n);
% stack=uint16(0);
% stack(n,n,planes,frames)=0;
% photons=photons/planes;
% n=ones(1,n);
for frame=1:frames
%    photon=photons(:,frame(n));
%    for plane=1:planes
%       stack(:,:,plane,frame)=poissrnd(x{plane}*(y{plane}.*photon));
%    end
   stack(:,:,:,frame) = imfilter(poissrnd(stack(:,:,:,frame)+b),psf,'symmetric')+gnoise*abs(randn(sy,sx,sz));
   fig=statusbar(frame/frames,fig);
%    if isempty(fig);
%       stack=stack(:,:,:,1:frame);
%       break;
%    end
end
delete(fig);

end

%One-dimensional PSF.
%
% function psf=gaussian(fwhm,centers,points)
% psf=-4*log(2)./fwhm.^2;
% if numel(psf) > 1
%    [~,psf]=ndgrid(points,psf);
% end
% [points,centers]=ndgrid(points,centers);
% psf=exp(psf.*(points-centers).^2);
% psf=psf.*repmat(1./sum(psf,1),size(psf,1),1);


%Intensity trace of an emitter (photons per frame).
%
function photons=brightness(Ion,Ton,Toff,Tbl,frames)
cycle=Ton + Toff;
cycles=10 + ceil(frames/cycle);
times=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
times(1)=times(1) - rand*(sum(times(1:10)));
times=cumsum(times(:));
while times(end) < frames
   cycles=ceil(2*(frames - times(end))/cycle);
   cycles=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
   cycles(1)=cycles(1) + times(end);
   times=[times;cumsum(cycles(:))];
end
times=times.';
Ton=times(2:2:end) - times(1:2:end);
Tbl=cumsum(Ton) + Tbl*log(rand);
n=find(Tbl > 0);
if any(n)
   Ton(n(2:end))=0;
   n=n(1);
   Ton(n)=Ton(n) - Tbl(n);
   times(2*n)=times(2*n) - Tbl(n);
end
photons=[zeros(size(Ton));Ion*Ton];
photons=cumsum(photons(:));
photons=diff(interp1(times,photons,0:frames,'linear',0));

end

% function photons=brightness(Ion,Ton,Toff,frames)
% cycle=Ton + Toff;
% cycles=10 + ceil(frames/cycle);
% times=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
% times(1)=times(1) - 5*cycle;
% times=cumsum(times(:));
% while times(end) < frames
%    cycles=ceil(2*(frames - times(end))/cycle);
%    cycles=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
%    cycles(1)=cycles(1) + times(end);
%    times=[times;cumsum(cycles(:))];
% end
% times=times.';
% photons=times(2:2:end) - times(1:2:end);
% photons=[zeros(size(photons));Ion*photons];
% photons=cumsum(photons(:));
% photons=diff(interp1(times,photons,0:frames));
