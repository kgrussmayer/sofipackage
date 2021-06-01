
% Select a results folder of main_MCSOFI2D
addpath('funcs')
addpath(['..',filesep,'..',filesep,'measureresolution',filesep,'funcs_ad']);

pname = uigetdir;

list = dir(pname);
list(1:2) = [];

for k = 1:length(list)
    if strfind(list(k).name,'_sofi_c1')
        ind = strfind(list(k).name,'_sofi_c1');
        tag = list(k).name(1:ind-1);
        ind = strfind(list(k).name,'.');
        fileExt = list(k).name(ind(end):end);
        break;
    end
end

sofi_c = [];
sofi_lin = [];

for k = 1:3
    sofi_c{k} = loadData([pname,filesep,tag,'_sofi_c',num2str(k),fileExt]);
end
for k = 1:3
    sofi_lin{k} = loadData([pname,filesep,tag,'_sofi_lin',num2str(k),fileExt]);
end


%% process and display the data

res_c = []; SNR_c = [];
for ch = 1:3
    for k = 1:size(sofi_c{ch},3)
        temp = sofi_c{ch}(:,:,k);
        temp = temp(1:min(size(temp,1),size(temp,2)),1:min(size(temp,1),size(temp,2)),:);
       [~,res,SNR] = getDcorr(apodImRect(temp,20));
       res_c(ch,k) = res;
       SNR_c(ch,k) = SNR;
    end
end

res_lin = []; SNR_lin = [];
for ch = 1:3
    for k = 1:size(sofi_lin{ch},3)
        temp = sofi_lin{ch}(:,:,k);
        temp = temp(1:min(size(temp,1),size(temp,2)),1:min(size(temp,1),size(temp,2)),:);
       [~,res,SNR] = getDcorr(apodImRect(temp,20));
       res_lin(ch,k) = res;
       SNR_lin(ch,k) = SNR;
    end
end

%% plot the results

figure(30000)
plot(res_c','linewidth',2); hold on
plot(res_lin','--','linewidth',2); hold off
title('Cut off frequency')
xlabel('Running parameter')
legend('raw c1', 'raw c2' , 'raw c3' , 'lin c1', 'lin c2' , 'lin c3')

figure(30001)
plot(SNR_c','linewidth',2); hold on
plot(SNR_lin','--','linewidth',2); hold off
title('SNR')
xlabel('Running parameter')
legend('raw c1', 'raw c2' , 'raw c3' , 'lin c1', 'lin c2' , 'lin c3')
%% save results

saveFigure(30000,pname,'resolution_estimate','tif')
saveFigure(30001,pname,'SNR_estimate','tif')

%%
figure(52345);
legends = [];
for k = 1:length(sofi_c)
    s =[];
    for h = 1:size(sofi_c{k},3)
        s(:,h) = getRadAvg(abs(fftshift(fftn(fftshift(sofi_c{k}(:,:,h))))));
    end
    kx = linspace(0,size(sofi_c{k},1)/size(sofi_c{1},1),size(s,1));
    plot(kx,log(mean(s,2)+1)); hold on
    legends{end+1} = ['Sofi ',num2str(k)];
end
hold off

legend(legends)
xlabel('Normalized frequencies')
ylabel('log of abs ')


k0 = 2*pi/0.6;
NA = 1.27;
dx = 0.108;
kmax = pi/dx;
ncut = NA*k0/kmax;

hold on; plot([ncut ncut],[0 5],'linewidth',2); hold off
