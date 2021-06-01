function val = clamp(val,minVal,maxVal)

if val < minVal
    val = minVal;
elseif val > maxVal 
       val = maxVal;
else
    return
end