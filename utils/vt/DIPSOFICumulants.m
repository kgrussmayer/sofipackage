function [ StackOUT,fwhms ] = DIPSOFICumulants( StackIN,S )



% EXTRACT DATA
Orders = S.Orders;
Nord = length(Orders);

% CALCULATE
[c,grid] = sofiCumulants2D(double(StackIN),1,[],[],Orders);
t = 1;
[c,fwhms] = sofiFlatten([],c,grid,Orders);

for ii = 1: Nord
    StackOUT{ii}(:,:,1) = c{ii};
end


end

