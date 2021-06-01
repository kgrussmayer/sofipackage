function s = setnestedfield(s, field, val)

fList = fieldnames(s);

for k = 1:length(fList)
    
    if ~isstruct(s.(fList{k}))
        if strcmp(field,fList{k})
           s.(field) = val;
           break
        end
    else
        subStruct = s.(fList{k});
        subfList = fieldnames(subStruct);
        for n = 1:length(subfList)
            if strcmp(field,subfList{n})
                subStruct.(field) = val;
                s.(fList{k}) = subStruct;
                break
            end
        end
    end
    
end