function showxsubs(desc,varargin)

figure('units','normalized','outerposition',[0 0 1 1],'visible','on'),
numPlots = nargin-1;

% choose layout of the subplots according to number of images
switch numPlots
    case 1
        rows=1;cols=1;
    case 2
        if isfield(desc(1),'vertical')
            if desc(1).vertical ==1
                rows=2;cols=1;
            else
                rows=1;cols=2;
            end
        else
            rows=1;cols=2;
        end
    case 3
        if isfield(desc(1),'vertical')
            if desc(1).vertical ==1
                rows=3;cols=1;
            else
                rows=1;cols=3;
            end
        else
            rows=1;cols=3;
        end
    case 4
        if isfield(desc(1),'vertical')
            if desc(1).vertical ==1
                rows=4;cols=1;
            else
                rows=2;cols=2;
            end
        else
            rows=2;cols=2;
        end
    case 5
        rows=2;cols=3;
    case 6
        rows=2;cols=3;
    otherwise
        rows = ceil(sqrt(numPlots));cols=rows;
end

% create all subplots with axis labels and appropriate colormaps
for jj=1:numPlots
    subplot(rows,cols,jj);imshow(varargin{jj},[]);
    try
        if isfield(desc(jj),'xlabel')
            xlabel(desc(jj).xlabel);
        end
        if isfield(desc(jj),'colormap')
            if strcmp(desc(jj).colormap,'none')
            else
                load([desc(jj).colormap,'.mat']);
                cmapName = desc(jj).colormap;
                colormap(eval(cmapName));freezeColors;
            end
        end
    catch
    end
    
end

% figure title
if isfield(desc(1),'title')
    title(desc(1).title);
end
% eof