[fn,pn]=uigetfile('*.tif','MultiSelect','on');
if ~iscell(fn)
    fn={fn};
end

stack=[];
for n=1:length(fn)
    fileinfo = imfinfo([pn fn{n}]);
    Nframes=length(fileinfo);

    for m=1:Nframes
        m
        if n*m==1
            stack=double(imread([pn fn{n}],m,'info',fileinfo));
        else
            stack(:,:,end+1)=double(imread([pn fn{n}],m,'info',fileinfo));
        end
    end
end

%%
stack=reshape(stack,size(stack,1),size(stack,2),Nframes,[]);
stack=permute(stack,[1 2 4 3]);

%% bin
M=stack;
mn=[size(M,1),size(M,2)];
binsize=2;
M=sum(reshape(M,binsize,[]));
M=permute(reshape(M,mn(1)/binsize,mn(2),[]),[2 1 3]);
M=sum(reshape(M,binsize,[]));
M=permute(reshape(M,mn(2)/binsize,mn(1)/binsize,[]),[2 1 3]);