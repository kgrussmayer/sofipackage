function [Rs] = radSum4(img,sector)
% Computes radial sum over a specified angle
%
% inputs: 
% in1 ...   input image matrix (2D matrix)
% sector ... angle for circular sum in rad, [start angle, end angle]
%            for example [0, pi/2]
% note: input image matrix img is assumed to be a 2D square matrix 
%
% ouput:
% Rs ... radial sum over an angle specified by [sector(1),sector(2)]
%
% tomas.lukes@epfl.ch


[N, M] = size(img);

dimMax = max(N,M);


%% Compute radial sum
% sector = [0,1*pi/6];
% sector = [1*pi/6,2*pi/6];
% sector = [2*pi/6,3*pi/6];
% sector = [3*pi/6,4*pi/6];
% sector = [4*pi/6,5*pi/6];
% sector = [5*pi/6,6*pi/6];
% sector = [6*pi/6,7*pi/6];
% sector = [7*pi/6,8*pi/6];
% sector = [9*pi/6,10*pi/6];
% sector = [10*pi/6,11*pi/6];
% sector = [11*pi/6,12*pi/6];

[X, Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);  % Cartesian grid
[theta, rho] = cart2pol(X, Y);                                  % Convert to polar coordinates
theta(1:find(theta(:,end)==0),:) = flip(abs(theta(1:find(theta(:,end)==0),:)),2)+pi;

rho = round(rho);
% theta = theta';
i = cell(floor(dimMax/2) + 1, 1);

if isempty(sector)
    for r = 1:floor(dimMax/2)
        i{r} = find(rho == r);
    end
else
    if length(sector)>2
        for r = 1:floor(dimMax/2)
            i{r} = find((rho == r & theta >= sector(1) & theta <= sector(2)) |...
                (rho == r & theta >= sector(3) & theta <= sector(4)));
        end
        
    else
        
        for r = 1:floor(dimMax/2)
            i{r} = find(rho == r & theta >= sector(1) & theta <= sector(2));
        end   
    end
end

Rs = zeros(1, floor(dimMax/2)+1);
% temp = zeros(size(img));

for r = 1:floor(dimMax/2)
%     temp(i{r}) = 1 ;
    Rs(1, r) = sum(img(i{r}));
end
% Rs = [img(rho==0),Rs]; % add zero frequency
% figure, imshow(temp)