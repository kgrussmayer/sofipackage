function id = getID(n)
if nargin < 1; n = 3;end

t = clock;

y = num2str(t(1));

m = num2str(t(2));
while length(m) < 2
    m = ['0',m];
end
d = num2str(t(3));
while length(d) < 2
    d = ['0',d];
end
h = num2str(t(4));
while length(h) < 2
    h = ['0',h];
end
mi = num2str(t(5));
while length(mi) < 2
    mi = ['0',mi];
end

switch n
    case 1
        id = y;
    case 2
        id = [y,m];
    case 3
        id = [y,m,d];
    case 4
        id = [y,m,d,'_',h];
    case 5
        id = [y,m,d,'_',h,mi];
    otherwise
        id = [y,m,d];
end
