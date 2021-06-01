function [imOut,relChangeAll] = deconvAugLag2D(im,h,settings)
% image deconvolution using the augmented Lagrangian approach with split-Bregman iterations

%
% Initialize
%
maxIter = settings.maxIter;
alpha = settings.alpha;
gamma = settings.gamma;
beta = settings.beta;
Lp = 1;% which Lp norm to use 
reltol = settings.reltol; % relative tolerace used as convergence criterion

[sy, sx, sz] = size(im);

%
% Prepare matrices
%

% size of the input image
imsize = size(im);

% FDx, FDy ... FFT of derivative operators along x and y 
% FDx = fft2([1 -1],imsize(1),imsize(2));
% FDy = fft2([1; -1],imsize(1),imsize(2));

FDx = repmat(fft2([1 -1],imsize(1),imsize(2)),[1 1 sz]);
FDy = repmat(fft2([1; -1],imsize(1),imsize(2)),[1 1 sz]);
if sz>1
%     d_z(1,1,1)= 1; d_z(1,1,2)= -1;
%     FDz = (fft2(d_z, imsize));
    FDz = repmat(fft2([1 -1],imsize(1),imsize(3)),[sx 1 1]);
    FDz = reshape(FDz,imsize(1),imsize(2),imsize(3));
else
    FDz = zeros(imsize);
end

% FH ... FFT of PSFs
% FH = psf2otf(h,imsize);
FH = h;

% FH=fft3d(im,h);
% FH = conj(fft2(hshift,imsize(1),imsize(2))).*fft2(h,imsize(1),imsize(2));
FHTH = conj(FH).*FH;

% FGs ... FFT of H^T*g
FGu = fftn(im);
FGs = conj(FH).*FGu;

DTD = conj(FDx).*FDx + conj(FDy).*FDy;
% DTD = conj(FDx).*FDx + conj(FDy).*FDy + conj(FDz).*FDz;
% FU ... FFT of u
FU = zeros(imsize);

% extra variables for Bregman iterations
Bx = zeros(imsize);
By = zeros(imsize);
% Bz = zeros(imsize);
Vx = zeros(imsize);
Vy = zeros(imsize);
% Vz = zeros(imsize);

%
% Iterative minimization of the cost function
%
for i = 1:maxIter
%     i
    FUp = FU; % store previous estimate
%     b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By)+conj(FDz).*fft2(Vz+Bz));     
    b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By)); 
    FU = b./(FHTH + beta/gamma*DTD);

    % Prepare the Lp prior
    Pr = asetupLnormPrior(Lp,alpha,beta);
    % get a new estimation of the auxiliary variable v
    xD = real(ifft2(FDx.*FU));
    yD = real(ifft2(FDy.*FU));
%     zD = real(ifft2(FDz.*FU));
    
    xDm = xD - Bx;
    yDm = yD - By;
%     zDm = yD - Bz;
    
    nDm = sqrt(sum(xDm.^2,3) + sum(yDm.^2,3));
    Vy = Pr.fh(yDm,nDm);
    Vx = Pr.fh(xDm,nDm);
%     Vx = beta*xDm/(2+beta);
%     Vy = beta*yDm/(2+beta);
%     Vz = beta*zDm/(2+beta);
    
    % update Bregman variables
    Bx = Bx + Vx - xD;
    By = By + Vy - yD;
%     Bz = Bz + Vz - zD;
    
    % evaluate the relative convergence criterion
    relChange = sqrt(sum(abs(FUp(:)-FU(:)).^2))/sqrt(sum(abs(FU(:)).^2));
    relChangeAll(i)=relChange;
    
    if relChange < reltol
        disp(['Stopped after ',num2str(i),' iterations. Convergence was reached within the predefined criteria.']);
        break;
    end
    
end
% disp(['Image deconvolved after ',num2str(i),' iterations. Final relative change is ', num2str(relChange)]);
imOut = abs(ifftn(FU));
    
end
