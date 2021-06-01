%[sofi,grid]=sofiCumulants3D(stack,first,frames,region,orders,gpu)
%-----------------------------------------------------------------
%
%Raw cross-cumulant images from image stacks.
%
%Inputs:
% stack     Image stacks [X Y Z frames]
% first     First image index                {1}
% frames    Number of images                 {all}
% region    Image region [pixel]             {{[1 width],[1 height]}}
% orders    Cross-cumulant orders            {1:4}
% gpu       Force or forbid CUDA usage       {auto}
%
%Outputs:
% sofi      Raw cross-cumulants
% grid      Partitions and pixels

%Copyright � 2012 Marcel Leutenegger et al, �cole Polytechnique F�d�rale de Lausanne,
%Laboratoire d'Optique Biom�dicale, BM 5.142, Station 17, 1015 Lausanne, Switzerland.
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
function [sofi,grid]=sofiCumulants2Dt(stacks,first,frames,region,orders,gpu)
[X,Y,Z,F]=size(stacks);
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
      x=int16(x(1):x(2));
      y=int16(y(1):y(2));
      roi=numel(x) < X || numel(y) < Y;
      X=numel(x);
      Y=numel(y);
   end
   if any([X Y frames] < 4)
      error('sofi:dimensions','The image stack is too small.');
   end
   if nargin < 5 || isempty(orders)
      orders=1:4;
   else
      orders=orders(0 < orders & orders < 7 & orders == fix(orders));
   end
   fig=statusbar(1/3,fig);
   grid=sofiGridsT(orders);           % 3D partitions and pixels
   disp('GridsT');
   tasks=(grid(1).pixels*[1;X;X*Y]).'; % corner pixel coordinates
   terms=grid(1).terms;                % terms pixel coordinates
   for T=1:numel(terms)
      terms{T}=[numel(terms{T}) sort(tasks(terms{T}))];
   end
   tasks=int32(cat(2,X*Y*Z,T,terms{:}));
   if roi
      stacks=stacks(x,y,:,first:first+frames-1);
   elseif frames < F
      stacks=stacks(:,:,:,first:first+frames-1);
   end
   means=sum(stacks,4)/frames;         % average for zero means
   fig=statusbar(2/3,fig);
   if nargin < 6 || isempty(gpu)
      gpu=cudaAvailable;
   end
   if gpu
      xpu=fileparts(mfilename('fullpath'));
      xpu=parallel.gpu.CUDAKernel(fullfile(xpu,'private','gpu.ptx'),fullfile(xpu,'private','gpu.cu'));
      %
      % If linking to the GPU kernel fails, try to recompile gpu.cu.
      % Change to the private folder and launch the Nvidia compiler.
      % The -arch option should be set to the compute capability of
      % your CUDA-capable graphics card. In the example, it is 2.1:
      %
      % >> !nvcc -arch=sm_21 -ptx gpu.cu
      %
      xpu.ThreadBlockSize=min(256,max(32*fix(sqrt(X*Y*Z)/32),32));
      xpu.GridSize=ceil(X*Y*Z/xpu.ThreadBlockSize(1));
      terms= gpuArray.zeros(X,Y,Z,T);
      tasks=gpuArray(tasks);
      means=gpuArray(means);
   else
      xpu=@cpu;
      terms=zeros(X,Y,Z,T);
   end
   fig=statusbar('Partial products...',fig);
   step=fix(frames/250);
   for frame=1:frames
      images=stacks(:,:,:,frame);
      if gpu
         images=gpuArray(images);
      end
      images=double(images);
      images=images - means;
%       gpuDevice(1)
      terms=feval(xpu,terms,images,tasks);
      if frame == frames || mod(frame,step) < 1
         fig=statusbar(frame/frames,fig);
         if isempty(fig)
            break;
         end
      end
   end
catch msg
   delete(fig);
   rethrow(msg);
end
if frame < frames
   error('sofi:abort','User abort resulted in non-zero means.');
else
   sofi=cumulants(grid,orders(orders > 1),frame,terms,fig,gpu);
   if any(orders == 1)
      sofi{1}=gather(means(2:end-2,2:end-2,:));
   end
end


%Check for CUDA-1.3-capable or newer graphics card.
%
function gpu=cudaAvailable
try
   gpu=parallel.gpu.GPUDevice.current();
   gpu=isa(gpu,'parallel.gpu.CUDADevice') && gpu.DeviceSupported;
catch msg
   gpu=false;
end


%Cross-cumulant images from partial product terms.
%
% grid      Partitions and pixels
% orders    Cumulant orders
% frames    Number of images
% terms     Partial products
% fig       Statusbar handle
%
% sofi      Cross-cumulant images
%
function sofi=cumulants(grid,orders,frames,terms,fig,gpu)
fig=statusbar('Cross cumulants...',fig);
fact=cumprod([1 -1:-1:1-max(orders(:))]/frames);
sofi=cell(size(grid));
s=size(terms);
s=[s(1:2)-3 s(3)];                     % first order size
x=cell(4,1);
y=cell(4,1);
z=cell(2,1);
for f=1:4
   x{f}=int16(f:f+s(1)-1);             % all shift indices
   y{f}=int16(f:f+s(2)-1);
end
z{1}=int16(1:s(3));
z{2}=int16([2:s(3) 1]);
if gpu
   fact=gpuArray(fact);
end
cur=0;
num=sum(cellfun(@numel,{grid(orders).parts})) + numel(orders);
for order=orders(:).'
   pixels={0};                         % cumulant groups
   pixels=pixels(ones(order,order,order));
   P=int16(1:numel(pixels));
   parts=grid(order).parts;
   shifts=grid(order).shifts;
   for p=1:numel(parts)
      dx=x(shifts{p}(:,:,1));          % shift indices
      dy=y(shifts{p}(:,:,2));
      dz=z(shifts{p}(:,:,3));
      part=parts{p};                   % term indices
      n=size(part,2);
      F=int16(1:n);
      for m=P
         term=fact(n);                 % partition
         for n=F
            term=term.*terms(dx{m,n},dy{m,n},dz{m,n},part(m,n));
         end
         pixels{m}=pixels{m} + term;   % cumulant
      end
      cur=cur + 1;
      fig=statusbar(cur/num,fig);
      if isempty(fig)
         return;
      end
   end
   term=zeros(s*order);                % full image
   for k=1:order
   for m=1:order
      for n=1:order
         term(k:order:end,m:order:end,n:order:end)=gather(pixels{k,m,n});
      end
   end
   end
   sofi{order}=term(1:1+end-order,1:1+end-order,1:1+end-order);
   cur=cur + 1;
   fig=statusbar(cur/num,fig);
   if isempty(fig)
      return;
   end
end
delete(fig);
