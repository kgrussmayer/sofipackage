% Author : Adrien Descloux 
% Use this function to plot the output of the SOFI processing

function plotSOFIstruct(c,figID)

if nargin < 2
    f = figure;
    figID = f;
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
set(f,'WindowScrollWheelFcn',@figScroll)

him = imagesc(c{1}(:,:,1,1)); % create an imagesc object
colorbar;

Ns = length(c);
N = size(c{1});

if length(N) == 4 
    Ny = N(1); Nx = N(2); Nc = N(3); Nw = N(4);
elseif length(N) == 3
    Ny = N(1); Nx = N(2); Nc = 1; Nw = N(3);
elseif length(N) == 2
    Ny = N(1); Nx = N(2); Nc = 1; Nw = 1;
end

set(gca,'position',[0.05 0.15 0.5 0.7])
if Nw > 1
    hs_w = uicontrol('style','slider','String','subw',...
            'units','normalized','position',[0.7 0.15 0.05 0.7],...
            'Callback',{@hs_w_Callback});
            uicontrol('style','text','String','Sub Window',...
            'units','normalized','position',[0.67 0.85 0.1 0.05]);
    set(hs_w,'Min',1,'max',Nw,'value',1,'SliderStep',[1/(Nw-1) 1/(Nw-1)])
    addlistener(hs_w,'Value','PreSet',@hs_w_Callback);
end
if Nc > 1
    hs_c = uicontrol('style','slider','String','channel',...
            'units','normalized','position',[0.8 0.15 0.05 0.7],...
            'Callback',{@hs_c_Callback});
    set(hs_c,'Min',1,'max',Nc,'value',1,'SliderStep',[1/(Nc-1) 1/(Nc-1)])
        uicontrol('style','text','String','Channel',...
            'units','normalized','position',[0.77 0.85 0.1 0.05]);
end

hs_o = uicontrol('style','slider','String','order',...
            'units','normalized','position',[0.9 0.15 0.05 0.7],...
            'Callback',{@hs_o_Callback});
set(hs_o,'Min',1,'max',Ns,'value',1,'SliderStep',[1/(Ns-1) 1/(Ns-1)])
        uicontrol('style','text','String','SOFI order',...
            'units','normalized','position',[0.87 0.85 0.1 0.05]);

h.ht = uicontrol('style','text','String',['SOFI order : ',num2str(1),', chanel : ',num2str(1),', sub window : ',num2str(1)],...
            'units','normalized','position',[0.05 0.85 0.5 0.05],'fontSize',12);

if length(size(c{1})) == 3
    set(hs_c,'visible','off')
    Nc = 1;
end

uicontrol('style','text','String',['Data size : # orders : ',num2str(length(c)),', # chanels : ',...
            num2str(Nc),', sub windows : ',num2str(Nw)],...
            'units','normalized','position',[0.05 0.93 0.6 0.05],'FontSize',13);

h.hb_m = uicontrol('style','pushbutton','String','Compute SOFI mean',...
                    'units','normalized','position',[0.2,0.03,0.2,0.05],...
                    'Callback',@computeMean);
h.hb_s = uicontrol('style','pushbutton','String','Compute STD per subsequence',...
                    'units','normalized','position',[0.6,0.03,0.2,0.05],...
                    'Callback',@computeSTD);
                
h.figID = figID;
h.him = him;
if Nw > 1
    h.hs_w = hs_w;
end
if Nc > 1
    h.hs_c = hs_c;
end
h.hs_o = hs_o;
h.c = c;
h.Nw = Nw;

guidata(figID,h);
update;
end

% MAIN UPDATE FUNCTION
function update
h = guidata(gcf);
order = round(get(h.hs_o,'value'));

N = size(h.c{order});
w = 1;
c = 1;
if length(N) == 4
    h.Nc = N(3);
    if h.Nc > 1
        set(h.hs_c,'Min',1,'max',h.Nc,'SliderStep',[1/(h.Nc-1) 1/(h.Nc-1)])
        c = round(get(h.hs_c,'value'));
        if c > size(h.c{order},3)
            c = size(h.c{order},3);
            set(h.hs_c,'value',c);
        end
    end
    if h.Nw > 1
        w = round(get(h.hs_w,'value'));
    end

    set(h.him,'CData',h.c{order}(:,:,c,w));
elseif length(N) == 3
%     h.Nc = 2;
    if h.Nc > 1
        set(h.hs_c,'Min',1,'max',1,'value',1,'SliderStep',[1 1])
        c = round(get(h.hs_c,'value'));
    end
    if h.Nw > 1
        w = round(get(h.hs_w,'value'));
    end
    set(h.him,'CData',h.c{order}(:,:,w));
elseif length(N) == 2
    set(h.him,'CData',h.c{order});
end
set(gca,'xLim',[1 N(2)],'yLim',[1 N(1)])

set(h.ht,'String',['SOFI order : ',num2str(order),', chanel : ',num2str(c),', sub window : ',num2str(w)])
guidata(h.figID,h)
end


function computeMean(source,event)
h = guidata(gcf);
order = round(get(h.hs_o,'value'));
N = size(h.c{order});
if length(N) == 4
    Nc = N(3);
elseif length(N) == 3
    Nc = 1;
elseif length(N) == 2
    Nc = 1;
end
prompt = {'Starting subwindow : ','End subwindow'};
name = 'Define mean range';
numlines = 1;
defaultanswer={'1',num2str(h.Nw)};
answer=inputdlg(prompt,name,numlines,defaultanswer);
figure
p = get(gcf,'position');
set(gcf,'position',round([p(1) p(2) p(3) p(4)]))
for k = 1:Nc
    if length(N) == 4
        subplot(ceil(Nc/ceil(sqrt(Nc))),ceil(sqrt(Nc)),k);imagesc(mean(h.c{order}(:,:,k,str2double(answer{1}):str2double(answer{2})),4));colorbar
        title(['SOFI order ',num2str(order),', Channel : ',num2str(k)])
    else
        subplot(ceil(Nc/ceil(sqrt(Nc))),ceil(sqrt(Nc)),k);imagesc(mean(h.c{order}(:,:,str2double(answer{1}):str2double(answer{2})),3));colorbar
        title(['SOFI order ',num2str(order)])
    end
    
end

set(gcf,'name',['Subwindows : ',answer{1},' to ',answer{2}])
end

function computeSTD(source,event)
h = guidata(gcf);
order = round(get(h.hs_o,'value'));
N = size(h.c{order});
s = [];
if length(N) == 4
    figure; temp = [];
   for j = 1:N(3)
       for k = 1:N(4)
        temp(k,j) = std(std(h.c{order}(:,:,j,k)));
       end
       s{j} = ['Sofi channel ',num2str(j)];
   end
   plot(temp);xlabel('Subsequence')
   legend(s)
else
    
end

end

function figScroll(source,event)
    dir = event.VerticalScrollCount;
    h = guidata(gcf);
    val = get(h.hs_w,'value');
    val = val-dir;
    
    if val > get(h.hs_w,'Max')
        val = get(h.hs_w,'Max');
    elseif val < 1
        val = 1;
    end
    set(h.hs_w,'value',val);
    update;
    guidata(h.figID,h)
end

function hs_w_Callback(source,event)
        update;      
end
function hs_c_Callback(source,event)
        update;      
end
function hs_o_Callback(source,event)
        update;      
end