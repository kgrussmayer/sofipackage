fname = 'E:\tlukes\data_HOF_JHI_sample\data_all\test';
% fname = [fname,filesep,'PALM_160513_pXJ41_CD4_wt_PSCFP2_Jurkat_fix_ser01_1x256_PLL001.tif-resultsdriftcor-80nm','.json'];
fname = [fname,filesep,'PALM_160513_pXJ41_CD4_wt_PSCFP2_Jurkat_fix_ser01_1x256_PLL001.tif-resultsdriftcor-105b','.json'];

data = loadjson(fname);

polyy = [];
polyx = [];

for ii = 1:data.yFunction.n
    disp(ii)
%     polyy = [polyy,polyval(data.yFunction.polynomials{ii}.coefficients,...
%         linspace(1,2,length(data.yFunction.knots(ii):data.yFunction.knots(ii+1))-1)-1)];
    polyy = [polyy,polyval(data.yFunction.polynomials{ii}.coefficients,...
        ones(1,data.yFunction.knots(ii+1)-data.yFunction.knots(ii)))];
    
    polyx = [polyx,polyval(data.xFunction.polynomials{ii}.coefficients,...
        ones(1,data.xFunction.knots(ii+1)-data.xFunction.knots(ii)))];
    
%     polyy = [polyy,polyval(data.yFunction.polynomials{ii}.coefficients,data.yFunction.knots(ii):data.yFunction.knots(ii+1)-1)];
end

polyy = [polyy,polyval(data.yFunction.polynomials{end}.coefficients,1)];
polyx = [polyx,polyval(data.xFunction.polynomials{end}.coefficients,1)]; 
% polyy = [polyy,polyval(data.yFunction.polynomials{end}.coefficients,data.yFunction.knots(ii+1))];

%%
figure, 
subplot(121);
plot(data.driftDataFrame,data.driftDataY,'g');hold on;
plot(1:data.maxFrame,polyy);ylim([-100 100]);
title(fname);

subplot(122);
plot(data.driftDataFrame,data.driftDataX,'r');hold on;
plot(1:data.maxFrame,polyx);ylim([-100 100]);

% plot(data.driftDataFrame,data.driftDataX)
