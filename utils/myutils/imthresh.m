function [out,threshold]=imthresh(im)
%IMTHRESH Iterative (optimal) image thresholding.
% CMP Vision Algorithms http://visionbook.felk.cvut.cz
%  
% Threshold a grayscale image using an automatic iterative threshold selection
% . 
%  
% Usage: [out,threshold] = imthresh(im)
% Inputs:
%  im   [m x n]  Input grayscale image with intensities 0... 255.
% Outputs:
%   out  [m x n]  Binary output image, contains ones where the pixel
%     values of im are above or equal to the threshold and
%     zeros elsewhere.
%  threshold  1x1  Threshold value.
% 



% We pre-calculate the histogram and cumulative histogram. This way
% we can evaluate in constant time the mean of the pixels above or below
% threshold. 
histogram = hist( im(:), 0:255 );
hist_times_gray = cumsum( histogram.*[0:255] );
cumulative_histogram = cumsum( histogram );

% A first approximation of the background mean mean_1 is  
% the mean of the corner pixels. A first approximation of the foreground
% mean mean_2 is the mean of the rest of the image. 
% The initial threshold
% is set as an average of the two means. 
[m,n] = size(im);
sum_background = sum( im([1 m  n*(m-1)+1  n*m]) );
num_pix_background = 4;
mean_1 = sum_background / num_pix_background;
mean_2 = (sum(im(:))-sum_background) / (n*m-num_pix_background);
threshold = ceil( (mean_1+mean_2)/2 );

if (threshold~=0)&&(cumulative_histogram(threshold)==0)
    threshold_old=threshold;
end 

% We calculate new foreground and background means defined by the 
% threshold and use them to calculate the new threshold, repeating until
% convergence. The procedure is very fast because the actual thresholding
% only takes place once the final threshold is found.
threshold_old = 0; 
while threshold~=threshold_old
  threshold_old = threshold;
  mean_1 = hist_times_gray(threshold) / cumulative_histogram(threshold);
  %WTF end means here
  mean_2 = (hist_times_gray(end)-hist_times_gray(threshold+1)) / ...
           (cumulative_histogram(end)-cumulative_histogram(threshold+1));
  threshold = ceil( (mean_1+mean_2)/2 );
end
out = im>=threshold;

