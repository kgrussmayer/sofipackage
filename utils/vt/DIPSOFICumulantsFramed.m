function [ StackOUT,fwhms] = DIPSOFICumulantsFramed( Stack,S )
% ORIGINAL FRAMED CALCULATION OF SOFI

[Y,X,Z] = size(Stack);
Orders = S.Orders;
Nord = length(Orders);
WSize = S.WSize;

NT = floor(Z/WSize);
fwhms = zeros(1,NT);

if(NT<1)
    warning(horzcat('Not enough samples: ',num2str(Z)));
end
for ns = 1:NT
    % CALCULATE
    SubStack = double(Stack(:,:,(ns-1)*WSize+1:(ns*WSize)));
    [c,grid] = sofiCumulants2D(SubStack,1,[],[],Orders);
    [c,fwhms(ns)] = sofiFlatten([],c,grid,Orders);

    for ii = 1:Nord
        StackOUT{ii}(:,:,ns) = c{ii};
    end
end



end

