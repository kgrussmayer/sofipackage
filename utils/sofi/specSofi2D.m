function mat_new = specSofi2D(mat)

% mat = S2_decn_seq;

intstep = 2; % interpolation step
mat_new = zeros(intstep*[size(mat,1),size(mat,2)]);

% pixels of the initial grid
% mat_c = abs(fft(mat,[],3)).*conj(abs(fft(mat,[],3)));
% mat_c = mean(mat_c(:,:,2:end),3);
mat_c = cpsd(mat,mat); 

indicesx = 1:intstep:size(mat_new,2);
indicesy = 1:intstep:size(mat_new,1);
mat_new(indicesy,indicesx)=mat_c;

%%% columns
% every second column a)
mat_c1 = mat(:,1:2:end,:);
mat_c2 = mat(:,2:2:end,:);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesx = 2:intstep*2:size(mat_new,2);
indicesy = 1:intstep:size(mat_new,1);
mat_new(indicesy,indicesx)=mat_c12;

% %every second column b)
mat_c1 = circshift(mat(:,1:2:end,:),-1,2);
mat_c2 = mat(:,2:2:end,:);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesx = 4:intstep*2:size(mat_new,2);
indicesy = 1:intstep:size(mat_new,1);
mat_new(indicesy,indicesx)=mat_c12;

% %%% rows
% every second row a)
mat_c1 = mat(1:2:end,:,:);
mat_c2 = mat(2:2:end,:,:);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesy = 2:intstep*2:size(mat_new,1);
indicesx = 1:intstep:size(mat_new,2);
mat_new(indicesy,indicesx)=mat_c12;

%every second row b)
mat_c1 = circshift(mat(1:2:end,:,:),-1,1);
mat_c2 = mat(2:2:end,:,:);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesy = 4:intstep*2:size(mat_new,1);
indicesx = 1:intstep:size(mat_new,2);
mat_new(indicesy,indicesx)=mat_c12;

%%% diagonals
% every second diagonal a)
mat_c1 = mat(1:2:end,1:2:end,:);
mat_c2 = mat(2:2:end,2:2:end,:);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesy = 2:intstep*2:size(mat_new,1);
indicesx = 2:intstep*2:size(mat_new,2);
mat_new(indicesy,indicesx)=mat_c12;

%every second diagonal b) - nekde je chyba
mat_c1 = circshift(mat(1:2:end,1:2:end,:),-1,1);
mat_c1 = circshift(mat_c1,-1,2);
mat_c2 = mat(2:2:end,2:2:end,:);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesy = 4:intstep*2:size(mat_new,1);
indicesx = 4:intstep*2:size(mat_new,2);
mat_new(indicesy,indicesx)=mat_c12;

% every second diagonal c)
mat_c1 = mat(1:2:end,2:2:end,:);
mat_c2 = circshift(mat(2:2:end,1:2:end,:),-1,2);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesy = 2:intstep*2:size(mat_new,1);
indicesx = 4:intstep*2:size(mat_new,2);
mat_new(indicesy,indicesx)=mat_c12;

%every second diagonal d)
mat_c1 = mat(2:2:end,1:2:end,:);
mat_c2 = circshift(mat(1:2:end,2:2:end,:),-1,1);
mat_c12 = cpsd(mat_c1,mat_c2); 

indicesy = 4:intstep*2:size(mat_new,1);
indicesx = 2:intstep*2:size(mat_new,2);
mat_new(indicesy,indicesx)=mat_c12;

function mat_c12 = cpsd(mat_c1,mat_c2)
fp = size(mat_c1,3);
% mat_c12 = abs(ifft(abs(fft(mat_c1,[],3)).*conj(abs(fft(mat_c2,[],3))),[],3));
mat_c12 = abs(fft(mat_c1,[],3)).*conj(abs(fft(mat_c2,[],3)));
% mat_c12 = mean(mat_c12(:,:,2:floor(end/2)),3); % integrate only over first half of samples
mat_c12 = mean(mat_c12(:,:,2:floor(end/4)),3); 
% mat_c12 = mean((real(ifft(mat_c12(:,:,1),[],3))),3); 
% mat_c12 = (real(ifft(mat_c12(:,:,1),[],3))); 

% mat_c12 = circshift(fftshift(real(ifft(mat_c12,2*fp,3)),3),-1);
% mat_c12 = circshift(real(ifft(mat_c12,2*fp,3)),-1);
% mat_c12 = mean(abs(mat_c12(:,:,fp:end-1)),3);

