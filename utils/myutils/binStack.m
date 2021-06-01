function stack2 = binStack(stack, blockx)
% blockx = 5;
stack2 = (zeros(floor(size(stack,1)./blockx),floor(size(stack,2)./blockx),size(stack,3)));

% % subtract background
% stack = stack -median(stack(:));

for ii = 1:blockx
    for jj = 1:blockx
        stackdec = stack(ii:blockx:end,jj:blockx:end,:);
        stackdec = stackdec(1:size(stack2,1),1:size(stack2,2),:);
        stack2 = stack2 + stackdec;
    end
end