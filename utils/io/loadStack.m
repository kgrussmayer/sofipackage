function [stack,frames] = loadStack(settings, str)
% load standard tif or big tif file

    stack = load_tifFile([settings.io.imageFile, str],settings.sys.sub,settings.io.roi);
    
    frames=size(stack,3);
    
    try
        if  size(stack,3) == 1 % the stack needs to be loaded as a big tiff file
            [stack,frames] = load_bigtif([settings.io.imageFile],settings.sys.sub);
        end 
    catch
            disp('load_bigtif failed');
    end
    
    if settings.io.ro ==1 % reorder the stack - for Nikon in Leuven
    stack2 = stack;    
        stack2(:,:,1:3:end)=stack(:,:,1:frames/3);  
        stack2(:,:,2:3:end)=stack(:,:,frames/3+1:2*frames/3);
        stack2(:,:,3:3:end)=stack(:,:,2*frames/3+1:frames);   
    stack = stack2;
    clear stack2;
    end