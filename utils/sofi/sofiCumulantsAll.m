function [sofi_c, settings,stats ] = sofiCumulantsAll(stack, settings)
% Calculate cumulants, perform flattening

    
    
    if settings.sys.frames > size(stack,3)
        settings.sys.frames = size(stack,3);
    end
    if settings.sys.wsize > size(stack,3)
        settings.sys.wsize = size(stack,3);
    end
    
    frames = size(stack,3);
    settings.sys.frames = frames;
    
    count=1;
    if settings.sys.jk ==0
        stats=0;
        for jj = 1:settings.sys.wsize:ceil(frames/settings.sys.wsize)*settings.sys.wsize
            if (jj+settings.sys.wsize-1) > size(stack,3) 
                substack = stack(:,:,jj:end);
                if size(substack,3) < 100
                    break;
                end
            else
                substack = stack(:,:,jj:jj+settings.sys.wsize-1);
            end
            tic
                if cudaAvailable==1
                [c,grid]=sofiCumulants2D(substack,1,[],[],settings.sys.orders); % cross cumulants, zero time lag 
            %     [cf,grid]=fourierCumulants(substack,1,[],[],orders);
        %         ck{2} = specSofi2D(substack);
        %         ck{2} = ksofi(substack);
                else
                [c,grid]=sofiCumulants(substack,1,[],[],settings.sys.orders);    
                toc
                end
%             end

            c=sofiAllFlatten(c,settings.sys.orders);
    %         [c,fwhm]=sofiFlatten([],c,grid,settings.sys.orders);
    %         fwhms(count)=fwhm;
            for io=settings.sys.orders
                sofi_c{io}(:,:,count)=medfilt2(c{io},[2 2]);
%                 sofi_c{io}(:,:,count)=c{io};
    %             sofi_ck{io}(:,:,count)=ck{io};
        %         sofi_cf{io}(:,:,count)=cf{io};
            end
            count = count+1;
        end
    else
        %@Vytautas: THIS LINE WAS WRONG, No-one checked the code thiw jk
        %estimator ON
        %for jj = settings.sys.start:settings.sys.wsize:ceil((frames-settings.sys.start+1)/settings.sys.wsize)*settings.sys.wsize
        for jj = 1:settings.sys.wsize:ceil(frames/settings.sys.wsize)*settings.sys.wsize
            if (jj+settings.sys.wsize-1) > size(stack,3)
                substack = stack(:,:,jj:end);
            else
                substack = stack(:,:,jj:jj+settings.sys.wsize-1);
            end
        if size(substack,3)>99
            tic
            [c,grid,bia,snr,var]=sofiCumulantsSNR(substack,1,[],[],settings.jk.orders); % cross cumulants, zero time lag 
    %         [c,grid,bia,snr,var]=sofiCumulantsSNRseq5_live(substack,1,[],[],settings.sys.orders,[],settings.sys.block);
            %     [cf,grid]=fourierCumulants(substack,1,[],[],orders);
    %         ck{2} = specSofi2D(substack);
    %         ck{2} = ksofi(substack);
            toc
         end

        c=sofiAllFlatten(c,settings.jk.orders);
%         [c,fwhm]=sofiFlatten([],c,grid,settings.sys.orders);
%         fwhms(count)=fwhm;
        for io=settings.jk.orders
            sofi_c{io}(:,:,count)=medfilt2(c{io},[2 2]);
            bias{io}(:,:,count)=bia{io};
            snrs{io}(:,:,count)=snr{io};
            vars{io}(:,:,count)=var{io};
%             sofi_c{io}(:,:,count)=c{io};
%             sofi_ck{io}(:,:,count)=ck{io};
    %         sofi_cf{io}(:,:,count)=cf{io};
        end
        count = count+1;
        end

  stats.b = bias;
  stats.s = snrs;
  stats.v = vars;
    end