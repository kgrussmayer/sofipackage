function out = isnestedfield(s, field)

fList = fieldnames(s);
out = 0;
for k = 1:length(fList)
    
    if ~isstruct(s.(fList{k}))
        if strcmp(field,fList{k})
            out = 1;
            break
        end
    else
        subStruct = s.(fList{k});
        subfList = fieldnames(subStruct);
        for n = 1:length(subfList)
            if strcmp(field,subfList{n})
                out = 1;
                break
            end
        end
    end
    
end