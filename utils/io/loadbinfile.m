function stack = loadbinfile(pn,fn,N)
% load stack of images from binary file, function assumes that binary
% images are in the same folder as "setting.mat" which contains information
% about acquisition
% Tomas Lukes, tomas.lukes@epfl.ch

fid1=fopen([pn fn]);
info=dir([pn fn]);

load([pn 'settings.mat']);
nx=cam1.ROIPosition(3);
ny=cam1.ROIPosition(4);

if nargin<3
    N=info.bytes/2/nx/ny;
end

k=0;
stack=uint16(zeros(ny,nx,N));
while k<N
    k=k+1;
    stack(:,:,k)=fread(fid1,[ny,nx],'*uint16');
    disp(['loading.. ',num2str(k)]);
end
fclose(fid1);
% eof