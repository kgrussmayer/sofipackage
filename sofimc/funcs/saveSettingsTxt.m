function saveSettingsTxt(settings,fname)

fList = fieldnames(settings);
fid = fopen( fname, 'wt' );

for k = 1:length(fList)
    
    if ~isstruct(getfield(settings,fList{k}))
        t = getfield(settings,fList{k});
        if isstr(t)
            fprintf(fid,'%s',['settings.',fList{k},' = ',t]);
        else
            fprintf(fid,'%s',['settings.',fList{k},' = ',num2str(t)]);
        end
        
            fprintf(fid,'\n');
    else
        fprintf(fid,'\n');
        subStruct = getfield(settings,fList{k});
        subfList = fieldnames(subStruct);
            for n = 1:length(subfList)
                t = getfield(subStruct,subfList{n});
                	if isstr(t)
                        fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ',t]);
                    else
                        fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ',num2str(t)]);
                    end
                    fprintf(fid,'\n');
            end
    end
    
end

fclose(fid);