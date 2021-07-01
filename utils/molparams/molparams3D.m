% Calculate molparams for 3D SOFI
% Tomas Lukes, tomas.lukes@epfl.ch

% settings
settings.molpar.thresh = 0.2;
sz = settings.sys.nplanes;

% calc molparames for each pair of consecutive z planes

for ii = 1:sz
    disp(ii);
    sofi2D{2} = sofizt{2}(:,:,1+2*(ii-1));
    sofi2D{3} = sofizt{3}(:,:,1+3*(ii-1));
    sofi2D{4} = sofizt{4}(:,:,1+4*(ii-1));
%     mparams = molparams(sofizt2D,sofi_lin,settings);
    
    sofilin2D{2} = sofid{2}(:,:,1+2*(ii-1));
    sofilin2D{3} = sofid{3}(:,:,1+3*(ii-1));
    sofilin2D{4} = sofid{4}(:,:,1+4*(ii-1));
    
    [ratio,density,brightness]=sofiParametersMy2(sofi2D);
    sofimask = sofilin2D{4};
    alpha = mean(sofimask,3);
    alpha_max = max(alpha(:));
    alpha = imadjust(alpha./alpha_max,[settings.molpar.thresh 1],[0 1]);

    % alpha = alpha./max(alpha(:));
    alpha(alpha>0) = 1;
    density(isinf(density)) = 0;

    bsofi=sofiBalance(sofilin2D,ratio);
    alpha = bsofi./max(bsofi(:));
    
    I_ratio = ratio.*alpha;
    I_bright = brightness.*alpha;
    I_density = density.*alpha;
    
    mparams = makestruct(ratio,brightness, density, I_ratio,I_bright,I_density, alpha, bsofi);
    
end
