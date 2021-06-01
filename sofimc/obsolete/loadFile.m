%%% load sofi files 
clear all;
filepath = 'C:\Users\Lukestom\Documents\ÈVUT\Ph.D. III.semestr\SOFI method\Original\SOFI_Measurements\007_316Hz\';

file1='data1.bin';
file2='data2.bin';

load([filepath 'settings.mat']);
% load([filepath '_myfile.mat']);
%%
fid1=fopen([filepath file1]);
fid2=fopen([filepath file2]);
info=dir([filepath file1]);

nx=cam1.ROIPosition(3);
ny=cam1.ROIPosition(4);

N=info.bytes/2/nx/ny;

k=0;
data1c=uint16(zeros(ny,nx,N));
data2c=data1c;
while k<N
    k=k+1;
    data1c(:,:,k)=fread(fid1,[ny,nx],'*uint16');
    data2c(:,:,k)=fread(fid2,[ny,nx],'*uint16');
end
clear k;
fclose(fid1);
fclose(fid2);
%%

data1c = data1c(:,1:236,:);
data2c = data2c(:,1:236,:);

%%
figure, imshow(mean(data1c,3),[])
figure, imshow(mean(data2c,3),[])