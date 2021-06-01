function saveSettingsTxt(settings,fname)

fList = fieldnames(settings);
fid = fopen( fname, 'wt' );

for k = 1:length(fList)
    
    if ~isstruct(settings.(fList{k}))
        t = settings.(fList{k});
        if isstr(t)
            fprintf(fid,'%s',['settings.',fList{k},' = ',t]);
        elseif iscell(t)
            fprintf(fid,'%s',['settings.',fList{k},' = ']);
            if numel(t) > length(t) || length(t) > 15
                fprintf(fid,'%s','TOO LARGE TO BE DISPLAYED');
            else
                for h = 1:length(t)
                    if isstr(t{1})
                        fprintf(fid,'%s',[t{h},', ']);
                    else
                        fprintf(fid,'%s',[num2str(t{h}),', ']);
                    end
                end
            end
        end
        
            fprintf(fid,'\n');
    else
        fprintf(fid,'\n');
        subStruct = settings.(fList{k});
        subfList = fieldnames(subStruct);
            for n = 1:length(subfList)
                t = subStruct.(subfList{n});
                	if isstr(t)
                        fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ',t]);
                    elseif iscell(t)
                        fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ']);
                        if numel(t) > length(t) || length(t) > 15
                            fprintf(fid,'%s','TOO LARGE TO BE DISPLAYED');
                        else
                            for h = 1:length(t)
                                if isstr(t{1})
                                    fprintf(fid,'%s',[t{h},', ']);
                                else
                                    fprintf(fid,'%s',[num2str(t{h}),', ']);
                                end
                            end
                        end
                    else
                        if size(t,1) == 1 || size(t,2) == 1
                            fprintf(fid,'%s',['settings.',fList{k},'.',subfList{n},' = ',num2str(t)]);
                        else
                            if size(t,1) >= size(t,2)
                                for h = 1:size(t,2)
                                    fprintf(fid,'%s\n',['settings.',fList{k},'.',subfList{n},' = ',num2str(t(:,h))]);
                                end
                            else
                              	for h = 1:size(t,1)
                                    fprintf(fid,'%s\n',['settings.',fList{k},'.',subfList{n},' = ',num2str(t(h,:))]);
                                end
                            end
                        end
                    end
                    fprintf(fid,'\n');
            end
    end
    
end

fclose(fid);