function descfile_sofi3D(settings)
% generate a text file with description of the data
% Tomas Lukes, tomas.lukes@epfl.ch

% Create content of the text file

fnames = fieldnames(settings);
clear A;
A{1} = sprintf('%s','Description file for data generated using 2D SOFI'); 
A{2} = sprintf('%s','----------------------------------------------------------------');
A{3} = sprintf('%s','');

ii = numel(A);ii = ii+1;

for ff = 1:numel(fnames)

    A{ii} = sprintf('%s',['Settings.',fnames{ff}]); ii = ii+1;
    subfnames = fieldnames(eval(['settings.',fnames{ff}]));
    for gg = 1:numel(subfnames)
        try
            A{ii} = sprintf('%s',['Settings.',fnames{ff},'.',subfnames{gg},' = ',...
                    num2str(getfield(settings,fnames{ff},subfnames{gg}))]);ii = ii+1; 
        catch err
        end
    end
    A{ii} = sprintf('%s','.....');ii = ii+1;
end
A{ii+1} = -1;
A = A';

% Write cell A into txt
outputFile = [settings.io.outputpath,filesep,'config.txt'];
% outputFile = ['config.txt'];

fid = fopen(outputFile, 'w');
for i = 1:ii
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\r\n', A{i});
    end
end
fclose(fid);

