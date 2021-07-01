function [sofi,grid]=mcSofiCumulantsSingle(stack,first,frames,region,orders,gpu)
%Calculate Raw cross cumulants images from image stack.
%
% Inputs:
% stack     Image stack
% first     First image index             {1}
% frames    Number of images              {all}
% region    Image region [pixel]          {{[1 width],[1 height]}}
% orders    Cross-cumulant orders         {1:4}
% gpu       Force or forbid CUDA usage    {auto}
%
% Outputs:
% sofi      Raw cross-cumulants
% grid      Partitions and pixels

% Copyright © 2014-2018 Marcel Leutenegger, Tomas Lukes
% École Polytechnique Fédérale de Lausanne,
% Laboratoire d'Optique Biomédical and
% Laboratory of Nanoscale Biology, http://lben.epfl.ch/
% tomas.lukes@epfl.ch

% This file is part of multicolorSOFI.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

[X,Y,F]=size(stack);

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
        orders=orders(0 < orders & orders < 9 & orders == fix(orders));
    end
    fig=statusbar(1/3,fig);
    grid=mcSofiGrids(orders,2);       % 2D partitions and pixels
    tasks=(grid(1).pixels*[1;X]).';     % corner pixel coordinates
    terms=grid(1).terms;                % terms pixel coordinates
    for T=1:numel(terms)
        terms{T}=[numel(terms{T}) sort(tasks(terms{T}))];
    end
    tasks=int32(cat(2,X*Y,T,terms{:}));
    
    if roi
        stack=stack(x,y,first:first+frames-1);
    elseif frames < F
        stack=stack(:,:,first:first+frames-1);
    end
    avg=sum(stack(:,:,:),3)/frames;  % average for zero means
    
    fig=statusbar(2/3,fig);
    if nargin < 6 || isempty(gpu)
        gpu=cudaAvailable;
    end
    if gpu
        xpu=fileparts(fileparts(mfilename('fullpath')));
        xpu=parallel.gpu.CUDAKernel(fullfile(xpu,'private','gpu.ptx'),fullfile(xpu,'private','gpu.cu'));
        %
        % If linking to the GPU kernel fails, try to recompile gpu.cu.
        % Change to the private folder and launch the Nvidia compiler.
        % The -arch option should be set to the compute capability of
        % your CUDA-capable graphics card. In the example, it is 2.1:
        %
        % >> !nvcc -arch=sm_21 -ptx gpu.cu
        %
        xpu.ThreadBlockSize=min(256,max(32*fix(sqrt(X*Y)/32),32));
        xpu.GridSize=ceil(X*Y/xpu.ThreadBlockSize(1));
        terms= gpuArray.zeros(X,Y,T);
        tasks=gpuArray(tasks);
        avg=gpuArray(avg);
    else
        xpu=@cpu;
        terms=zeros(X,Y,T);
    end
    fig=statusbar('Partial products...',fig);
    step=fix(frames/250);
    for frame=1:frames
        img=stack(:,:,frame);
        if gpu
            img=gpuArray(img);
        end
        img=double(img);
        img=img - avg;
        terms=feval(xpu,terms,img,tasks);
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

sofi=cumulants(grid,orders(orders > 1),frame,terms,fig,gpu);
if any(orders == 1)
    sofi{1}=gather(avg(2:end-2,2:end-2));
end


%Check for CUDA-1.3-capable or newer graphics card.
function gpu=cudaAvailable
try
    gpu=parallel.gpu.GPUDevice.current();
    gpu=isa(gpu,'parallel.gpu.CUDADevice') && gpu.DeviceSupported;
catch
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
s=s(1:2) - 3;                          % first order size
x=cell(4,1);
y=cell(4,1);
for f=1:4
    x{f}=int16(f:f+s(1)-1);             % all shift indices
    y{f}=int16(f:f+s(2)-1);
end
if gpu
    fact=gpuArray(fact);
end
cur=0;
num=sum(cellfun(@numel,{grid(orders).parts})) + numel(orders);
for order=orders(:).'
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
                term=term.*terms(dx{m,n},dy{m,n},part(m,n));
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
    for m=1:order
        for n=1:order
            term(m:order:end,n:order:end)=gather(pixels{m,n});
        end
    end
    sofi{order}=term(1:1+end-order,1:1+end-order);
    cur=cur + 1;
    fig=statusbar(cur/num,fig);
    if isempty(fig)
        return;
    end
end
delete(fig);