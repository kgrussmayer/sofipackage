function [stack,frames] = loadStack(settings, str)
% load standard tif or big tif file

[stack, frames] = load_tiff([settings.io.imageFile, str], ...
    settings.sys.sub, settings.io.roi);

if settings.io.ro ==1 % reorder the stack - for Nikon Leuven setup
    stack2 = stack;
    stack2(:,:,1:3:end)=stack(:,:,1:frames/3);
    stack2(:,:,2:3:end)=stack(:,:,frames/3+1:2*frames/3);
    stack2(:,:,3:3:end)=stack(:,:,2*frames/3+1:frames);
    stack = stack2;
    clear stack2;
end
% eof