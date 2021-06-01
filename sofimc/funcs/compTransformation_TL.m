function compTransformation_TL(hObject, eventdata, handles)
fig=get(hObject,'parent');
ax=findobj(fig,'type','axes');
axe1=ax(2);
axe2=ax(1);
cp1=findobj(axe1,'type','line','Color',[1 1 0]);
cp2=findobj(axe2,'type','line','Color',[1 1 0]);
im1=findobj(axe1,'type','image');
im1=get(im1,'CData');
im2=findobj(axe2,'type','image');
im2=get(im2,'CData');
xcp1=get(cp1,'XData');
ycp1=get(cp1,'YData');
xcp2=get(cp2,'XData');
ycp2=get(cp2,'YData');

% compute transformation map
pts1 = [xcp1',ycp1'];
pts2 = [xcp2',ycp2'];
assignin('base','pts1',pts1);
assignin('base','pts2',pts2);

tform=cp2tform([xcp1',ycp1'],[xcp2',ycp2'],'similarity');
assignin('base','tform',tform);

im1t=imtransform(im1,tform,'XData', [1 size(im2,2)], 'YData', [1 size(im2,1)], ...
                'Size', size(im2));

figure;imshowpair(sqrt(double(im1t)),sqrt(double(im2)),'blend');axis xy;axis equal;axis tight;

end
