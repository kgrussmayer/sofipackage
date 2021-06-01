function [stack,results] = bleachcor_monoexp(stack,settings,results)
%  perform simple bleaching correction based on monoexponential fiting
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings
% output:   stack ... bleaching corrected sequence of images
     
if ~isfield(results.blcor,'fitcurve')
    % mono exponential fitting algorithm 
    results = plotmitrace(stack,results);
end
mitrace_cor = zeros(size(results.blcor.mitrace)); % bleaching corrected meantrace
frames = size(stack,3);
for ii = 1:frames
    stack(:,:,ii) = stack(:,:,ii)/(results.blcor.fitcurve(ii));
    mitrace_cor(ii) = results.blcor.mitrace(ii)/(results.blcor.fitcurve(ii));
end
results.blcor.mitrace_cor = mitrace_cor;
% [corela] = stackcorel(stack,settings.blcor.MaxCorrSamp);
% results.blcor.corela_cor = corela;

    
