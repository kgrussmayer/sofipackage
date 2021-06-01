function [gcor,rho,coefs,parmap] = kest_v6(sofic,num,img1,stack)
% estimate rho_on and gamma correction from second and third order
% Tomas Lukes, tomas.lukes@epfl.ch
% sofic = sofi_c;
% img1 = wfbcg;

orders = 2:numel(sofic);
[~,~,k]=size(sofic{2});

for k = 1:k
    img2 = sofic{2}(:,:,k); 
    img3 = sofic{3}(:,:,k);
    [X,Y,~]=size(img2);
    img2=(img2(1:2:X,1:2:Y,:)); % take only subgrid corresponding to physical pixels
    
    [X,Y,~]=size(img3);
    img3=(img3(1:3:X,1:3:Y,:)); 

    % choose the optimal points
    [imbw1] = getbw(img1);
    [imbw2] = getbw(img2);
    [imbw3] = getbw(img3);
    
    mask = imbw1 + imbw2 + imbw3;
    ind = find(mask==3);
%     if num > length(ind)/2; 
%         num = round(length(ind)/2);
%         disp(['num: ',num2str(num)]);
%     end;
    
    mask = max(0,mask-ones(size(mask))*2);
    optpts = img2.*mask; % matrix of candidates for optimal points
    
    % discard points outside the range (0.3 - 0.9)*Ipeak
%     optpts(optpts < 0.3*max(optpts(:)) | optpts > 0.9*max(optpts(:))) = 0;
        temp = img2; 
%         temp = optpts;
%         temp(temp< 0.2*max(temp(:)) | temp > 0.8*max(temp(:))) = 0; 
%         figure, imshow(optpts,[])
        cb = 10;
        fwhm = 4;
        temp=LaplacianOfGaussian(temp,fwhm);
        temp([1:cb,end-cb],[1:cb,end-cb]) = 0;
       
    for ii = 1:num % take num biggest values
%         disp(ii)
        [~,ind] = max(temp(:));
        [my,mx]=ind2sub(size(img2),ind);

%         [~,ind] = max(optpts(:));
%         [my,mx]=ind2sub(size(optpts),ind);
        r = 3;
%         block_val = temp(my-r:my+r,mx-r:mx+r);
        [block_x,block_y] = meshgrid(-r:+r,-r:+r);
%         my-r:my+r,mx-r:mx+r
        block_d = sqrt(block_x.^2+block_y.^2);
        block_d(block_d<=r) = 1;
        block_d(block_d>r) = 0;
        
%         i1 = block_d.*img1(my-r:my+r,mx-r:mx+r); %g1
%         i2 = block_d.*img2(my-r:my+r,mx-r:mx+r); %g2
%         i3 = block_d.*img3(my-r:my+r,mx-r:mx+r); %g3
%         i1 = mean(i1(:));
%         i2 = mean(i2(:));
%         i3 = mean(i3(:));
        i1 = img1(my,mx); %g1
        i2 = img2(my,mx); %g2
        i3 = img3(my,mx); %g3

%         rho1(ii) = ((1/3^(3/2))*i2.^2 - (1/2^3)*i3.*i1)./...
%                     (2*(1/3^(3/2))*i2.^2 - (1/2^3)*i3.*i1);
%                 
%         rho1(ii) = abs(rho1(ii));
        temp(my,mx)=0;
%         
%         N = (1-rho1(ii))*i1.*i1./(2^(2/3)*rho1(ii)*i2);
%         N(isnan(N)) = 0;
%         parmap(ii,k).N = N;
        parmap(ii,k).mx = mx;
        parmap(ii,k).my = my;
        
%         parmap(ii,k).rho1 = rho1(ii);

        trace = squeeze(stack(my,mx,:));
%         trace = squeeze(mean(mean(repmat(block_d,1,1,size(stack,3)).*stack(my-r:my+r,mx-r:mx+r,:),1),2));
        
%         frames = length(trace);
        maxlag = 500;
        trace = xcorr(trace,trace,maxlag,'unbiased');%./(length(trace)*150);
        trace = trace(maxlag:end)./(maxlag);
%         [~,imin] = min(trace);
%         trace = trace(1:100)
        frames = length(trace);
        
        s = fitoptions('Method','NonlinearLeastSquares',...
                       'Lower',[0,0,0],...
                       'Upper',[Inf,Inf,Inf],...
                       'Startpoint',[1 1 1]);
        f = fittype('a*exp(-x/b)+c','options',s);
        
        [c2,~] = fit((1:frames-1)',trace(2:end),f);
        a = c2.a;
        b = c2.b;
        c = c2.c;
        
        % Show the fitting
%         figure;x=(0:frames-1);plot(x',trace);hold on;plot(x,a*exp(-x/b)+c,'k');
        
%         s = fitoptions('Method','NonlinearLeastSquares',...
%                        'Lower',[0,0,0],...
%                        'Upper',[Inf,Inf,Inf],...
%                        'Startpoint',[1 1 1]);
%         f = fittype('a*exp(-x/b)+c','options',s);
%         
%         [c2,~] = fit((0:frames-1)',trace,f);
%         a = c2.a;
%         b = c2.b;
%         c = c2.c;
%         
%         % Show the fitting
%         x=(0:frames-1);
%         exfit = a*exp(-x/b)+c;
%         figure;plot(x',trace);hold on;plot(x,exfit,'k');
        
        
        N = sqrt(a); 
        parmap(ii,k).N = N;
        
        a = i2.^(3/2);
        b = i3;
        c = 2^(9/4)/3^(3/2);
%         c = c*1/sqrt(N);
        ck = 1;
            
        c = c*1/(sqrt(ck)*sqrt(N));
            
        rho1b(ii) =(-sqrt(4*a*a*b*b*c*c+b^4) + 4*a*a*c*c + b*b)/(2*(4*a*a*c*c+b*b)); 
        rho2b(ii) =(sqrt(4*a*a*b*b*c*c+b^4) + 4*a*a*c*c + b*b)/(2*(4*a*a*c*c+b*b)); 
        parmap(ii,k).rho1b = rho1b(ii);
        parmap(ii,k).rho2b = rho2b(ii);
    end
%     rho1(rho1<0) = [];
    rho1 = mean(rho1b);
    rho(k) = rho1;
    
    for io = orders
        imgx = sofic{io}(:,:,k); 
        [X,Y,~]=size(imgx);
        imgx=(imgx(1:io:X,1:io:Y,:)); % take only subgrid corresponding to physical pixels
%         imgx(imgx< 0.2*max(imgx(:)) | imgx > 0.8*max(imgx(:))) = 0;
        
        for ii = 1:num
        [~,ind] = max(imgx(:));
        [my,mx]=ind2sub(size(imgx),ind);

        i1 = imgx(my,mx); %g2

        switch io
            case 2
                temp = abs(rho(k)*(1-rho(k)));
            case 3
                temp = abs(rho(k)*(1-rho(k))*(1-2*rho(k)));
            case 4
                temp = abs(rho(k)*(1-rho(k))*(1-6*rho(k)+6*rho(k)*rho(k)));
            case 5
                temp = abs(rho(k)*(1-rho(k))*(1-2*rho(k))*(12*rho(k)*rho(k)-12*rho(k)+1));
            case 6 
                temp = abs(rho(k)*(1-rho(k))*(120*rho(k)^4-240*rho(k)^3+150*rho(k)^2-30*rho(k)+1));
            otherwise 
                disp('Supported order has to be in the range [2,6]');
        end
        coefs(io,k) = temp;      
        gcor(io,k,ii) = log10(i1/temp)/log10(i1); % correction for the io order
%         gcor(io,k,ii) = log10(i1/(temp*N*(ck/io^(3/2))))/log10(i1); % correction for the io order
        imgx(my,mx)=0;
        end  
    end
end
gcor = mean(gcor,3);