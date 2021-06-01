function [z,CV] = SmoothPCA(y, Lambda)

% *************************************************************************
% Whittaker smoother for one-dimensional signals
%
% Input:
%   y:      One-dimensional signal (sampled at equal intervals)
%   Lambda: Smoothing parameter
%
% Output:
%   z:      Smoothed signal
%   CV:   	Leave-one-out cross-validation
%
% Based on: 
%   P. H. C. Eilers (2003). A Perfect Smoother. Analytical Chemistry, 
%   Vol. 75, No. 14, pp. 3631--3636.
%
% *************************************************************************

% Smoothing of the signal
m = length(y);
E = speye(m);
D = diff(speye(m), 2);
z = (E + Lambda * D' * D) \ y;

% Leave-one-out cross-validation
if nargout > 1
    
    % Exact hat diagonal
	H = inv(E + Lambda * D' * D);
    h = diag(H);
    
    % Cross-validation value
    r = (y - z) ./ (1 - h);
    CV = sqrt(r' * r / m);
    
end
