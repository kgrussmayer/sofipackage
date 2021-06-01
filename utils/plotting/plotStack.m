% plotStack(data,figID)
% ---------------------------------------
%
% Opens a basic gui that allows simplified navigation through 3D data stacks
% if figID is specified, first close(figID), then create figure(figID)
% if figID is not specified, open a new figure
%
% Inputs:
%   data        Data stack [X Y Z]
%   figID       Figure ID,  if set, close(figID),then create figure(figID)
%                           if not set, open new figure
%
% Outputs:
%   s           Structure containing experimental and processing parameters
%
% ---------------------------------------
% A detailled description of the theory supporting this program can be found in : 
% "Descloux, A., et al. "Combined multi-plane phase retrieval and 
%  super-resolution optical fluctuation imaging for 4D cell microscopy." 
%  Nature Photonics 12.3 (2018): 165."
%
%   Copyright © 2018 Adrien Descloux - adrien.descloux@epfl.ch, 
%   École Polytechnique Fédérale de Lausanne, LBEN/LOB,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.

function plotStack(data,figID)

dim = 2;

[s1,s2,s3] = size(data);
dat_min = min(data(:)); dat_max = max(data(:));
if nargin == 1
f = figure('visible','off','position',[360 500 400 285],'WindowScrollWheelFcn',@figScroll);
else
try close(figID);end
f = figure(figID);
set(f,'visible','off','position',[360 500 400 285],'WindowScrollWheelFcn',@figScroll);
end


set(f,'name',inputname(1))
hslide = uicontrol('style','slider','String','z',...
            'position',[330 30 30 200],...
            'Callback',{@slider_Callback});
set(hslide,'Min',1);set(hslide,'Max',s3);
set(hslide,'value',round((s3)/2));set(hslide,'SliderStep',[1/(s3) 1/(s3)]);
addlistener(hslide,'Value','PreSet',@slider_Callback);

htxt = uicontrol('style','edit','string','0',...
    'position',[330 235 30 15],'Callback',@edit_Callback);
hcmax = uicontrol('style','edit','string',num2str(dat_max),...
    'position',[335 250 50 8],'Callback',@hcmax_Callback,'enable','off');
hcmin = uicontrol('style','edit','string',num2str(dat_min),...
    'position',[335 260 50 8],'Callback',@hcmax_Callback,'enable','off');
hcmapbox = uicontrol('style','checkbox','value',1,'string','Colomap auto',...
                'Callback',@cmap_Callback,'position',[335 268 50 15]);
            
hpopcmap = uicontrol('style','popupmenu','string',{'parula','jet','hot','winter',...
                    'gray','copper','morgenstemning','isolum'},'units','normalized','position',[.1 .85 0.1 0.1],...
                    'value',1,'Callback',@popcmap_Callback);
            
hcbox = uicontrol('style','checkbox','value',0,'string','Scale colormap to 3D data',...
                'Callback',@cbocx_Callback,'position',[200 260 80 20]);

hpop = uicontrol('style','popupmenu','string',{'X','Y','Z'},'value',2,...
                'Callback',@pop_Callback,'position',[280 245 30 30]);

                 
ha = axes('Units','pixels','Position',[50,60,260,185]);

% by default, plot the midle of the data
switch dim
    case 1
        him = imagesc(reshape(data(round(end/2),:,:),[s2 s3]));
        
    case 2
        him = imagesc(reshape(data(:,round(end/2),:),[s1 s3])');
        title('Sliding along the "Y" direction')
    case 3
        him = imagesc(data(:,:,round(end/2)));
        colorbar
        title('Sliding along the "Z" direction');
    otherwise
        dim  = 3;
        disp('Display direction set to "z"')
        him = imagesc(data(:,:,round(end/2)));
        title('Sliding along the "Z" direction')
end

colormap('default')
colorbar
caxis('auto')

% define figure handles
handles.data = data;
handles.dim = dim;
handles.s = [s1 s2 s3];
handles.txt = htxt;
handles.slide = hslide;
handles.cbox = hcbox;
handles.cmapbox = hcmapbox;
handles.im = him;
handles.ax = ha;
handles.dmin = dat_min; handles.dmax = dat_max;
handles.hcmin = hcmin; handles.hcmax = hcmax;
handles.pop = hpop;
handles.popcmap = hpopcmap;
handles.f = f;

guidata(handles.f,handles)

% activation of the gui
set(f,'toolbar','figure');
set(f,'Units','normalized');
set(ha,'Units','normalized');
set(hslide,'Units','normalized');
set(hpop,'Units','normalized');
set(hcbox,'Units','normalized');
set(htxt,'Units','normalized');
set(hcmapbox,'Units','normalized');
set(hcmax,'Units','normalized');
set(hcmin,'Units','normalized');
set(f,'position',[0.0891    0.3417    0.3885    0.5037]);
set(f,'Visible','on');

end

function cbocx_Callback(source,~)
%     edit_Callback(source)
    slider_Callback(source);
end

function cmap_Callback(source,~)
h = guidata(gcf);
val = get(h.cmapbox,'value');
if val
    set(h.hcmax,'enable','off')
    set(h.hcmin,'enable','off')
else
	set(h.hcmax,'enable','on')
    set(h.hcmin,'enable','on')
end
slider_Callback(source);
end

function hcmax_Callback(source,~)
h = guidata(gcf);
cmax = str2double(get(h.hcmax,'string'));
cmin = str2double(get(h.hcmin,'string'));
cmax = clamp(cmax,h.dmin,h.dmax); cmin = clamp(cmin,h.dmin,h.dmax);
if cmax <= cmin; cmax = h.dmax;cmin = h.dmin; end
set(h.hcmax,'string',num2str(cmax))
set(h.hcmin,'string',num2str(cmin))
    slider_Callback(source);
end

function slider_Callback(source,~)

h = guidata(gcf);
h.dim = get(h.pop,'value');

% update slide properties
set(h.slide,'max',h.s(h.dim));
slid_val = get(h.slide,'value');
if slid_val > h.s(h.dim)
	set(h.slide,'value',h.s(h.dim))
end

val = get(h.slide,'value');
vmax = get(h.slide,'max');vmin = get(h.slide,'min');
val = round(val);
set(h.slide,'value',round(val));
% set(h.txt,'string',num2str(2*(round(val)-1)/(h.s(h.dim)-1)-1));
set(h.txt,'string',num2str(round(val)))

val = round((vmin + vmax)-val);
switch get(h.pop,'value')
    case 1
        temp = reshape(h.data(val,:,:),[h.s(2) h.s(3)])';
        set(h.im,'Cdata',temp)
        set(h.ax,'xlim',[0.5 h.s(2)+0.5],'ylim',[0.5 h.s(3)+0.5])
        title('Sliding along the "X" direction')
        xlabel('"Y" direction sdf');ylabel('"Z" direction')
    case 2
        temp = reshape(h.data(:,val,:),[h.s(1) h.s(3)])';
        set(h.im,'Cdata',temp)
        set(h.ax,'xlim',[0.5 h.s(1)+0.5],'ylim',[0.5 h.s(3)+0.5])
        title('Sliding along the "Y" direction')
        xlabel('"X" direction');ylabel('"Z" direction')
    case 3
        temp = h.data(:,:,val) ;
        set(h.im,'Cdata',temp)
        set(h.ax,'xlim',[0.5 h.s(2)+0.5],'ylim',[0.5 h.s(1)+0.5])
        title('Sliding along the "Z" direction')
        xlabel('"Y" direction');ylabel('"X" direction')
end

try
    cmapS = get(h.popcmap,'string');
    colormap(cmapS{get(h.popcmap,'value')})
end

if get(h.cbox,'value')
    if (h.dmin-h.dmax)
        caxis([h.dmin h.dmax])
    else
        caxis([h.dmin h.dmax+1])
    end
else
    if get(h.cmapbox,'value')
        caxis('auto')
    else
        cmax = str2double(get(h.hcmax,'string'));
        cmin = str2double(get(h.hcmin,'string'));
        caxis([cmin cmax])
    end
end

% save global changes
guidata(h.f,h);

end

function pop_Callback(source,event)
    h = guidata(gcf);
    set(h.txt,'string','0');
    h.dim = get(h.pop,'value');
    guidata(gcf,h);
    edit_Callback(source,event);
%     slider_Callback(source,event);
end

function popcmap_Callback(source,event)

    edit_Callback(source,event);
%     slider_Callback(source,event);
end

function edit_Callback(source,event)
    h = guidata(gcf);
    val = str2num(get(h.txt,'string'));
    if val < 1
%         val = floor(((val+1)/2)*(h.s(h.dim)-1)+1)
        val = floor(linmap(val,-1,1,1,h.s(h.dim)));
    end
        
    if val > h.s(h.dim)
    	val = h.s(h.dim);
    elseif val < 1
        val = 1;
    end
    
set(h.slide,'value',val);
slider_Callback(source,event);

end

function figScroll(source,event)

dir = event.VerticalScrollCount;

h = guidata(gcf);
val = get(h.slide,'value');

% val = round(val - dir.*0.05.*h.s(h.dim));
val = val -dir;
if val > h.s(h.dim)
    val = h.s(h.dim);
elseif val < 1
    val = 1;
end
set(h.slide,'value',val);
slider_Callback(source,event);

end

function  rsc = linmap(val,valMin,valMax,mapMin,mapMax)

% convert the input value between 0 and 1
tempVal = (val-valMin)./(valMax-valMin);

% clamp the value between 0 and 1
map0 = tempVal < 0;
map1 = tempVal > 1;
tempVal(map0) = 0;
tempVal(map1) = 1;

% rescale and return
rsc = tempVal.*(mapMax-mapMin) + mapMin;

end
