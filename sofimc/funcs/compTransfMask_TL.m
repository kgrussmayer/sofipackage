function compTransfMask_TL(im1,im2)
%%% transforms im1 and im2 images
%%% transformation matrix is stored in the work space in the tform variable
%%% lukestom 24.7. 2014

% segmentation
[y1,x1]=findpeaks(double(im1));
[y2,x2]=findpeaks(double(im2));

nx=size(im1,1);
ny=size(im1,2);
if nx>ny
    ysize=round(500*ny/nx);
    xsize=500;
else
    xsize=round(500*nx/ny);
    ysize=500;
end

%
s=get(0,'ScreenSize');
temp.mainFrameSize = [1400 800];
s=s(3:4)- temp.mainFrameSize;  %calculates the bottom left corner by leaving enough space on the top for the figure-header

fig=figure('Colormap',gray(256),'DoubleBuffer','on','Menubar','none','Name','compTranfMask','Tag','compTranfMask','NumberTitle','off','Position',[s/2 temp.mainFrameSize],'Resize','on');
axe1=axes('parent',fig,'units','pixels','position',[30 30 xsize ysize]);
axe2=axes('parent',fig,'units','pixels','position',[60+xsize 30 xsize ysize]);
go=uicontrol('parent',fig,'string','Calculate tranformation matrix','position',[30 670 300 60],'callback',@compTransformation_TL);
isource=imagesc(im1,'parent',axe1,'buttondownfcn',{@MousePressed});
itarget=imagesc(im2,'parent',axe2,'buttondownfcn',{@MousePressed});
sourcepeaks=line(x1,y1,'Linestyle','none','Marker','.','Color','red','Visible','off','Parent',axe1);
targetpeaks=line(x2,y2,'Linestyle','none','Marker','.','Color','red','Visible','off','Parent',axe2);
cpsource=line(1,1,'Linestyle','none','Marker','.','Color',[1 1 0],'Visible','on','Parent',axe1);
cptarget=line(1,1,'Linestyle','none','Marker','.','Color',[1 1 0],'Visible','on','Parent',axe2);
cpsn=text(1,1,'','Color',[1 1 0],'parent',axe1);
cptn=text(1,1,'','Color',[1 1 0],'parent',axe2);
set(cpsource,'xdata',[],'ydata',[]);
set(cptarget,'xdata',[],'ydata',[]);
disp('end');