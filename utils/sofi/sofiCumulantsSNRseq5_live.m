%[sofi,grid,bias,snrs,vars]=sofiCumulantsSNR(stack,first,frames,region,orders,gpu)
%---------------------------------------------------------------------------------
%
%Raw cross-cumulant images and jackknife statistics from image stack. If
%the user interrupts the calculation of the cross-cumulants, the jackknife
%statistics is calculated on the image frames digested before interruption.
%
%Inputs:
% stack     Image stack
%       or  TIFF file name
% first     First image index             {1}
% frames    Number of images              {all}
% region    Image region [pixel]          {{[1 width],[1 height]}}
% orders    Cross-cumulant orders         {1:4}
% gpu       Force or forbid CUDA usage    {auto}
%
%Outputs:
% sofi      Raw cross-cumulants
% grid      Partitions and pixels
% bias      Jackknife bias estimations
% snrs      Jackknife signal-to-noise ratios
% vars      Jackknife variance estimations

%Copyright © 2014 Marcel Leutenegger et al, École Polytechnique Fédérale de Lausanne,
%Laboratoire d'Optique Biomédicale, BM 5.142, Station 17, 1015 Lausanne, Switzerland.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [sofi,grid,bias,snrs,vars]=sofiCumulantsSNRseq5_live(stack,first,frames,region,orders,gpu,blockx)
file=ischar(stack);
if file
   tif=imfinfo(stack);
   X=tif(1).Height;                    % image dimensions
   Y=tif(1).Width;
   F=numel(tif);
   tif={'Info',tif};
else
   [X,Y,F]=size(stack);
end
try
   fig=statusbar('Initialization...');
   if nargin < 2 || isempty(first)
      first=1;
   end
   if nargin < 3 || isempty(frames)    % number of images
      frames=F + 1-first;
   end
   roi=nargin > 3 && ~isempty(region); % region defined?
   if roi
      x=max(1,min(region{1},X));
      y=max(1,min(region{2},Y));
      x=int16(min(x):max(x));
      y=int16(min(y):max(y));
      roi=numel(x) < X || numel(y) < Y;
      X=numel(x);
      Y=numel(y);
      if file
         tif{3}='PixelRegion';
         tif{4}={x([1 end]),y([1 end])};
      end
   end
   if any([X Y frames] < 4)
      error('sofi:dimensions','The image stack is too small.');
   end
   if nargin < 5 || isempty(orders)
      orders=1:4;
   end
   fig=statusbar(1/3,fig);
   grid=sofiGrids(orders,2,1);        % partitions and pixels
   orders=find(~cellfun(@isempty,{grid.dists}));
   if ~file
      if roi
         stack=stack(x,y,first:first+frames-1);
      elseif frames < F
         stack=stack(:,:,first:first+frames-1);
      end
   end
   fig=statusbar(2/3,fig);
   p=grid(1).pixels;                   % corner pixel coordinates
   m=cell(size(p));
   n=cell(size(p));
   F=int16(1:size(p,1));
   for f=F                             % indices limited to image size
      m{f}=int16(min(X,1+p(f,1):p(f,1)+X));
      n{f}=int16(min(Y,1+p(f,2):p(f,2)+Y));
   end
   t=grid(1).terms;                    % terms pixel coordinates
   P=int16(1:numel(t));
   facts=cell(size(m));                % factors storage
   if nargin < 6 || isempty(gpu)
      gpu=X*Y >= 65536 && cudaAvailable;
   end
   gpus={@gpu1 @gpu2 @gpu3 @gpu4 @gpu5 @gpu6 @gpu7 @gpu8};
   if gpu
      for f=F
         m{f}=gpuArray(int32(m{f}));   % use native data type
         n{f}=gpuArray(int32(n{f}));
      end
      terms={gpuArray(0)};
   else
      terms={0};
   end
   terms=terms(ones(size(t)));         % terms storage

   % Partial products for all frames.
   %
   fig=statusbar('Partial products...',fig);
   for frame=1:frames
      if file
         img=imread(stack,'Index',first+frame-1,tif{:});
      else
         img=stack(:,:,frame);
      end
      if gpu
         img=gpuArray(img);
      end
      img=double(img);
      for f=F
         facts{f}=img(m{f},n{f});      % all factors
      end
      if gpu
         for p=P
            terms{p}=arrayfun(gpus{numel(t{p})},terms{p},facts{t{p}});
         end
      else
         for p=P
            terms{p}=feval(gpus{numel(t{p})},terms{p},facts{t{p}});
         end
      end
      fig=statusbar(frame/frames,fig);
      if isempty(fig)
         break;
      end
   end
   frames=frame;                       % limit frames
   
   % Partial products and delete-1 cumulants.
   %
   fig=statusbar('Jackknife statistics...',fig);
   for p=P
      terms{p}=-terms{p};
   end
   termi=cell(size(terms));            % delete-1 terms

   sofi1=cell(size(grid));
   sofi1(orders)={0};                  % average
   sofi2=sofi1;                        % variance
   
   for frame=1:frames
      if file
         img=imread(stack,'Index',first+frame-1,tif{:});
      else
         img=stack(:,:,frame);
      end
      if gpu
         img=gpuArray(img);
      end
      img=double(img);
      for f=F
         facts{f}=img(m{f},n{f});      % all factors
      end
      if gpu
         for p=P
            termi{p}=arrayfun(gpus{numel(t{p})},terms{p},facts{t{p}});
         end
      else
         for p=P
            termi{p}=feval(gpus{numel(t{p})},terms{p},facts{t{p}});
         end
      end
      sofi=cumulants(grid,orders,frames-1,termi,gpu);
      
        for p=orders
        
        blockSize = blockx*p;

        temp = (binStack(abs(sofi{p}), blockSize));
        sofi1{p}= sofi1{p} + temp;
        sofi2{p}= sofi2{p} + temp.^2;
         
        end

      fig=statusbar(frame/frames,fig);
      if isempty(fig)
         break;
      end
   end
catch me
   delete(fig);
   rethrow(me);
end
sofi=cumulants(grid,orders,frames,terms,gpu);
sofisub = sofi;
for p=orders
    blockSize = blockx*p;
    sofisub{p} = (binStack(abs(sofi{p}), blockSize));
end
clear('term*');
delete(fig);

% Jackknife estimators.

if frame < frames
   warning('sofi:abort','User abort: jackknife estimation on %d of %d frames.',frame,frames);
end
bias=sofisub;                             % cumulants offsets
vars=sofisub;                             % cumulants variances
snrs=sofisub;                             % signal-to-noise ratios
for p=orders
   sofi1{p}=sofi1{p}/frame;
   sofi2{p}=sofi2{p}/frame;
   bias{p}=(frame-1)*(sofi1{p} - sofisub{p});
   vars{p}=(frame-1)*(sofi2{p} - sofi1{p}.^2);
   snrs{p}=abs(sofisub{p})./sqrt(vars{p});
end


%Check for CUDA-1.3-capable or newer graphics card.
%
function gpu=cudaAvailable
try
   gpu=parallel.gpu.GPUDevice.current();
   gpu=gpu.DeviceSupported;
catch
   gpu=false;
end


%Cross-cumulant images from negated partial product terms.
%
% grid      Partitions and pixels
% orders    Cumulant orders
% frames    Number of images
% terms     Partial products
% fig       Statusbar handle
%
% sofi      Cross-cumulant images
%
function sofi=cumulants(grid,orders,frames,terms,gpu)
fact=cumprod([-1 1:max(orders(:))-1]/frames);
sofi=cell(size(grid));
s=size(terms{1}) - 3;                  % first order size
x=cell(4,1);
y=cell(4,1);
for f=1:4
   x{f}=int16(f:f+s(1)-1);             % all shift indices
   y{f}=int16(f:f+s(2)-1);
end
if gpu
   fact=gpuArray(fact);
end
for order=orders
   pixels={0};
   pixels=pixels(ones(order));         % cumulant groups
   P=int16(1:numel(pixels));
   parts=grid(order).parts;
   shifts=grid(order).shifts;
   for p=1:numel(parts)
      dx=x(shifts{p}(:,:,1));          % shift indices
      dy=y(shifts{p}(:,:,2));
      part=parts{p};                   % term indices
      n=size(part,2);
      F=int16(1:n);
      for m=P
         term=fact(n);                 % partition
         for n=F
            term=term.*terms{part(m,n)}(dx{m,n},dy{m,n});
         end
         pixels{m}=pixels{m} + term;   % cumulant
      end
   end
   term=zeros(s*order);                % full image
   for m=1:order
      for n=1:order
         term(m:order:end,n:order:end)=gather(pixels{m,n});
      end
   end
   sofi{order}=term(1:1+end-order,1:1+end-order);
end
