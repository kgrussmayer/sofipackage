% Author : Adrien Descloux

function plotSOFI(sofi,figID)

if nargin < 2
    f = figure;
else
    try
        p = get(figID,'position');
        delete(figID)
        f = figure(figID);
        set(f,'position',p)
    catch
        f = figure(figID);
        set(f,'position',[100 100 800 500])
    end

end

Nc = length(sofi);
for k = 1:length(sofi)
    subplot(ceil(Nc/ceil(sqrt(Nc))),ceil(sqrt(Nc)),k);imagesc(sofi(k).imUp);
    colorbar
    title(sofi(k).name)
end


% grab zoom and pan handles
h.hZoom = zoom(gcf);
h.hPan = pan(gcf);

set(h.hZoom,'ActionPostCallback',{@updateZoom});
set(h.hPan,'ActionPostCallback',{@updateZoom});

guidata(gcf,h);
end

function updateZoom(obj, eventData)
ax = gca;
xl = ax.XLim; yl = ax.YLim;
f = gcf;
for k = 1:length(f.Children)
    set(f.Children(k),'XLim',xl);
    set(f.Children(k),'YLim',yl);
end

end
