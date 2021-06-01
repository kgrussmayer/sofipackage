function [outstruct] = makestruct(varargin)
% create a structure ("outstruct") which contains input variables, structure fields are
% given by the names of the input variables

outstruct = [];
for k = 1:length(varargin);
    varname = inputname(k);
    eval(['outstruct.',varname,'=','varargin{',num2str(k),'}']);
end

