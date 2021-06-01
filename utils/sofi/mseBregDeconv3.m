function [f] = mseBregDeconv3(im,h,maxIter,fcpx,lambda,nss,gamma)
% image deconvolution using Gaussian noise model and split-Bregman iterations
% im ... single or multi 2D or 3D input image
% maxIter ... maximum number of iterations
% fcpx ... cut-off frequency according to the optical system 
% related to the sampling frequency of the image - value in the range (0 - 1)
% lambda ... regularization parameter, will be estimated automaticly if it is not specified
% nss ... number of subsequences (for single 2D and single 3D nss = 1)
% gamma ... coefficient for linearization
% f ... output deconvolved, denoised image

%%% Tomas Lukes 
%%% last change 15.4.2015

%
% Initialize
%
if nargin < 6 || ~any(nss);
    nss = 1;
end
if nargin < 5 || ~any(lambda)
    lambda = estimateNoise(mean(mean(abs(im),4),3)/max(abs(im(:))));
    disp('Lambda parameter should be in the range (0.001 - 0.1)'); 
    disp(['Lambda estimate: ', num2str(lambda)]);
end
if lambda < 0.01; lambda = 0.001; end;
if lambda > 0.1; lambda = 0.1; end;

reltol = 1e-6; % relative tolerace used as convergence criterion
coeff_A = 0.1; % coeff_A > 0 => algoritm uses full norm of the regularization instead of seminorm

psf_size = size(h);

[sy,sx,sz] = getDimensions(im,nss);

masksign = ones(sy,sx,sz);
masksign (im<0)=-1;
im = abs(im);
if ~any(gamma)
im = (im./(max(im(:)))); 
else
im = (im./(max(im(:)))).^(gamma);   
end
im = masksign.*im; % preserve sign
    
im = bcgsub(im,sy,sx,sz); % subtract bacground

% pad borders
[im,Mask_im] = addBorders(im,sy,sx,sz,psf_size);

% get the dimensions of the padded image
[sy,sx,sz] = getDimensions(im,nss);

%
% Prepare matrices
%

% FH ... FFT of PSF
FH = psf2otf(h,[sy,sx,sz]);

% FDx, FDy ... FFT of derivative operators along x and y 
FDx = repmat(fft2([1 -1],sy,sx),[1 1 sz]);
FDy = repmat(fft2([1; -1],sy,sx),[1 1 sz]);

if sz>1  
    FDz = repmat(fft2([1 -1],sx,sz),[sy 1 1]);
    FDz = reshape(FDz,sy,sx,sz);
else
    FDz = zeros(sy,sx,sz);
end

% Fourier transform input images
dim = length(size(im)); 
if dim ==4
    % stack of 3D input images
    for ii = 1:nss
    FGu(:,:,:,ii) = fftn(im(:,:,:,ii));
    end
    FGu = sum(FGu,4);
        
elseif dim ==3 && nss == 1
    % 2D or 3D input image
    FGu = fftn(im);  
else
    % stack of 2D input images
    for ii = 1:nss
    FGu(:,:,ii) = fft2(im(:,:,ii));
    end 
    FGu = sum(FGu,3);
end

FHTH = conj(FH).*FH;

% FGs ... FFT of H^T*g
% FGs = conj(FH).*((abs(FGu).*spectralMask(sx,sy,sz)).*exp(1i*angle(FGu)));
FGs = conj(FH).*((abs(FGu)).*exp(1i*angle(FGu)));

DTD = conj(FDx).*FDx + conj(FDy).*FDy + conj(FDz).*FDz;

I = coeff_A*ones(sy,sx,sz);

% initial estimate for the f sub-problem
if dim ==2
    f = conj(FH).*FGu;
elseif dim==3 && nss ==1
    f = FH.*FGu; 
else
    f = FH.*fftn(mean(im,dim));
end
% f = zeros(sy,sx,sz); % iniate f with zeros

% extra variables for Bregman iterations
y1 = zeros(sy,sx,sz);
y2 = zeros(sy,sx,sz);
y3 = zeros(sy,sx,sz);

u1 = zeros(sy,sx,sz);
u2 = zeros(sy,sx,sz);
u3 = zeros(sy,sx,sz);

%
% main loop 
%

for i = 1:maxIter
    
    % f-subproblem
    f_old = f; % store previous estimate
    b = FGs + lambda*(conj(FDx).*fftn(u1-y1) + conj(FDy).*fftn(u2-y2)+conj(FDz).*fftn(u3-y3));     
    f = b./(I+ nss*FHTH + lambda*DTD);

    f = apodization(f,fcpx);
    f = real(ifftn(f));
    f = f./max(f(:));
   
    masksign = ones(sy,sx,sz);
    masksign (f<0)=-1;
    f = abs(f);
    
    f = imadjust(f,[lambda, 1],[0 1], 1); % suppress values below the uncertainty of the estimate
    
    f = masksign.*f;
    f = fftn(f); 

    % u-subproblem
    Dfx = Mask_im.*real(ifftn(FDx.*f));
    Dfy = Mask_im.*real(ifftn(FDy.*f));
    Dfz = Mask_im.*real(ifftn(FDz.*f));
    
    v1 = Dfx + y1;
    v2 = Dfy + y2;
    v3 = Dfz + y3;

    v  = sqrt(v1.^2 + v2.^2 + v3.^2);
   
    m=v;
    m(v==0) = 1; 
    v  = max(v - lambda, 0)./m; 
    
    u1 = v1.*v;
    u2 = v2.*v;
    u3 = v3.*v;
    
    % update Bregman variables
    y1 = y1 - u1 + Dfx;
    y2 = y2 - u2 + Dfy;
    y3 = y3 - u3 + Dfz;
       
    % evaluate the relative convergence criterion
    relChange = sqrt(sum(abs(f_old(:)-f(:)).^2))/sqrt(sum(abs(f(:)).^2));
    % estimate incremental SNR (iSNR)
    ISNR=calcSNR(real(ifftn(f_old)),real(ifftn(f)));
    
    disp(['Iteration:  ',num2str(i),' Relative change: ', num2str(relChange),'   ISNR: ',num2str(ISNR)]);
    
    if relChange < reltol
        disp(['Stopped after ',num2str(i),' iterations. Convergence was reached within the predefined criteria.']);
        break;
    end
    
end
disp(['Image deconvolved after ',num2str(i),' iterations. Final relative change is ', num2str(relChange)]);

f = apodization(f,fcpx);

f = real(ifftn(f));

f = imadjust(abs(f),[lambda, 1],[0 1], 1); 
% suppress values below the uncertainty of the estimate using `three-sigma rule' 
% for a zero-mean Gaussian random variable

% take out extra borders 
f = cutBorders(f,psf_size);

end

function [sy,sx,sz] = getDimensions(im,numIm)
    dim = length(size(im)); 

    if dim==4
        [sy,sx,sz,~] = size(im);
    elseif dim==3 && numIm ==1
        [sy,sx,sz] = size(im);
    else
        [sy,sx,~] = size(im);
        sz = 1;
    end

end

function [im, Mask_im] = addBorders(im,sy,sx,sz,psf_size)
    % padding the image according to PSF size (to suppress the border problem)

    Mask_im = ones(sy,sx,sz);
    Mask_im = padarray(Mask_im,psf_size);

    im = padarray(im,psf_size,'symmetric');
end

function f = cutBorders(f,psf_size)
    if size(f,3)>1 
    f = f(psf_size(2)+1:end-psf_size(2),psf_size(1)+1:end-psf_size(1),psf_size(3)+1:end-psf_size(3),:);
    else
    f = f(psf_size(2)+1:end-psf_size(2),psf_size(1)+1:end-psf_size(1));
    end
end

function im = apodization(im,omega)

    [sy,sx,sz] = size(im);
        
    if size(im,3)>1 
        [k_z, k_y]=meshgrid(-sz/2+1:sz/2,-sy/2+1:sy/2);
        k_r = sqrt(k_z.^2+k_y.^2);
        k_max = omega*max(k_r(:));
        bhs = cos(pi*k_r/(2*k_max));
        indi =  k_r > k_max ;
        bhs(indi) = 0;

        % apodize yz planes
        for ii = 1:size(im,2)
            yzplane = squeeze(im(:,ii,:));
            yzplane = fftshift(yzplane);

            % apply apodization
            yzplane = yzplane.*bhs; 
            im(:,ii,:)=ifftshift(yzplane);
        end

        [sy,sx,sz] = size(im);
        [k_z, k_x]=meshgrid(-sz/2+1:sz/2,-sx/2+1:sx/2);
        k_r = sqrt(k_z.^2+k_x.^2);
        k_max = omega*max(k_r(:));
        bhs = cos(pi*k_r/(2*k_max));
        indi =  k_r > k_max ;
        bhs(indi) = 0;

        % apodize xz planes
        for ii = 1:size(im,1)
            xzplane = squeeze(im(ii,:,:));
            xzplane = fftshift(xzplane);

            % apply apodization
            xzplane = xzplane.*bhs; 
            im(ii,:,:)=ifftshift(xzplane);
        end
    end
    
    [k_x, k_y]=meshgrid(-sx/2+1:sx/2,-sy/2+1:sy/2);
    k_r = sqrt(k_x.^2+k_y.^2);
    k_max = omega*max(k_r(:));
    bhs = cos(pi*k_r/(2*k_max));
    indi =  k_r > k_max ;
    bhs(indi) = 0;
    
    % apodize xy planes
    for ii = 1:size(im,3)
        xyplane = squeeze(im(:,:,ii));
        xyplane = fftshift(xyplane);

        % apply apodization
        xyplane = xyplane.*bhs; 
        im(:,:,ii)=ifftshift(xyplane);
    end
    
end

function im = bcgsub(im,sy,sx,sz)
% subtract bacground
    for ii = 1:sz 
        im(:,:,ii)=imtophat(im(:,:,ii),strel('disk',round(sqrt(sx^2+sy^2)/10)));
    end
end

function lambda=estimateNoise(im)

[sy, sx]=size(im);
im=double(im);

% compute sum of absolute values of Laplacian
M=[1 -2 1; -2 4 -2; 1 -2 1];
lambda=sum(sum(abs(conv2(im, M))));

% properly scale sigma
lambda=lambda*sqrt(0.5*pi)./(6*(sx-2)*(sy-2));
end

function m = spectralMask(sx,sy,sz)

fx = linspace(-1,1,sx);
fy = linspace(-1,1,sy);
[Fx,Fy] = meshgrid(fx,fy);

[THETA,RHO]=cart2pol(Fx,Fy); 

fc = 0.03;
m = RHO<=fc; 

m = double(m);

fwhm = fc*sx;
sigma = fwhm/(2*sqrt(2*log(2)));

g = fspecial('gaussian',[3*ceil(fwhm), 3*ceil(fwhm)],sigma);
m = imfilter(m,g);
m = m./max(m(:));
m = imcomplement(m);
m = repmat(m,[1 1 sz]);
end
