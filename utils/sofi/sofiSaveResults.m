function sofiSaveResults(settings,varargin)
% save SOFI images as tif stacks

avgsum = 'mean'; % calculate mean or sum to get the final image {mean, sum}

for order = settings.sys.orders

    % save all variables in varargin
    for k = 1:length(varargin);
        sofi = varargin{k};
        temp = feval(avgsum,sofi{order},3);
        varname = inputname(k+1); %first function variable is not in varargin
        saveImageStack(temp,settings.io.outputpath,[settings.io.imageName,varname,num2str(order)],[],16)
%       saveImageStack(sofi_c{order},settings.io.outputpath,[settings.io.imageName,'sofi_rawMovie',num2str(order)],[],16)
    end

end