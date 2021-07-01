function  rsc = linmap(val,valMin,valMax,mapMin,mapMax)

if nargin == 3
    mapMax = valMax;
    mapMin = valMin;
    valMin = min(val(:));
    valMax = max(val(:));
end

% convert the input value between 0 and 1
tempVal = (val-valMin)./(valMax-valMin);

% clamp the value between 0 and 1
map0 = tempVal < 0;
map1 = tempVal > 1;
tempVal(map0) = 0;
tempVal(map1) = 1;

% rescale and return
rsc = tempVal.*(mapMax-mapMin) + mapMin;
% eof