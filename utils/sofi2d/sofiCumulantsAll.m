function [sofi_c,settings,stats] = sofiCumulantsAll(stack,settings)
% Calculate cumulants, perform flattening

if settings.sys.frames > size(stack,3), settings.sys.frames = size(stack,3); end
if settings.sys.subseqlength > size(stack,3), settings.sys.subseqlength = size(stack,3); end
if isempty(settings.sys.start), settings.sys.start = 1; end
frames = settings.sys.frames;

count=1;
if settings.sys.jk == 0
    stats=0;
    for jj = settings.sys.start: ...
            settings.sys.subseqlength: ...
            ceil((settings.sys.frames-settings.sys.start+1)...
            /settings.sys.subseqlength)*settings.sys.subseqlength
        if (jj+settings.sys.subseqlength-1) > size(stack,3)
            substack = stack(:,:,jj:end);
        else
            substack = stack(:,:,jj:jj+settings.sys.subseqlength-1);
        end
        
        tic
        disp(size(substack))
        if cudaAvailable==1
            [c,grid]=sofiCumulants2D(substack,1,[],[],settings.sys.orders); % cross cumulants, zero time lag
        else
            [c,grid]=sofiCumulants(substack,1,[],[],settings.sys.orders);
        end
        toc
        
        c=sofiAllFlatten(c,settings.sys.orders);
        for io=settings.sys.orders
            sofi_c{io}(:,:,count)=medfilt2(c{io},[2 2]);
        end
        count = count+1;
    end
else
    for jj = settings.sys.start: ...
            settings.sys.subseqlength: ...
            ceil((settings.sys.frames-settings.sys.start+1)...
            /settings.sys.subseqlength)*settings.sys.subseqlength
        if (jj+settings.sys.subseqlength-1) > size(stack,3)
            substack = stack(:,:,jj:end);
        else
            substack = stack(:,:,jj:jj+settings.sys.subseqlength-1);
        end
        
        tic
        [c,grid,bia,snr,var]=sofiCumulantsSNR(substack,1,[],[],settings.jk.orders); % cross cumulants, zero time lag
        toc
        
        c=sofiAllFlatten(c,settings.jk.orders);
        for io=settings.jk.orders
            sofi_c{io}(:,:,count)=medfilt2(c{io},[2 2]);
            bias{io}(:,:,count)=bia{io};
            snrs{io}(:,:,count)=snr{io};
            vars{io}(:,:,count)=var{io};
        end
        count = count+1;
    end
    
    stats.b = bias;
    stats.s = snrs;
    stats.v = vars;
end
% eof