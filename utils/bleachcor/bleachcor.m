function [stack,results] = bleachcor(stack,settings,results)
%  wraper function for various bleaching correction algorithms 
%
% input:    stack ... sequence of input images acquired by camera
%           settings ... struct with all the processing settings
% output:   stack ... bleaching corrected sequence of images


switch settings.blcor.type
   
    case 'monoexp'
        % mono exponential fitting algorithm 
        [stack,results] = bleachcor_monoexp(stack,settings,results);
    
    case 'iir'
        % IIR filtering
        stack = bleachcor_iir(stack,settings);
        
    case 'PCA'
        % PCA based bleaching correction from Lille
        stack = bleachcor_pca(stack,settings);
    
end

%%
% figure, 
% imshow(mean(stack2,3),[]);