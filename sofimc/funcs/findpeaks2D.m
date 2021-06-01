function [segx,segy,segA,segE]=findpeaks2D(Id)

Is=LaplacianOfGaussian(Id,2.1);
% disp(['mean: ' num2str(mean(Is(:)))]);
% disp(['min: ' num2str(min(Is(:)))]);
% disp(['max: ' num2str(max(Is(:)))]);

%         figure;imagesc(Is);colorbar;
% disp(min(Is(:)));

%subtract background
% out.bgmap=double((Is<min(Is(:))*bgth));
%         figure;hist(Is(:),100);
bgmap=double((Is<min(Is(:))/10));
Is=bgmap.*Id;

% imagesp(Id,'Orig');
%         figure;imagesc(Is);colorbar;

%segment processed image
if any(Is(:))
    [S,n]=segmentImage(Is);
    % imagesp((S>0),'Segments');

    % imagesp(S,'Segments');

    %analyze segments
    [segx,segy,segA,segE]=analyzeSegments(S,Id,n,1); %ms=1 (Minimum subtraction on)
else
    segx = [];
    segy = [];
    segA = [];
    segE = [];
end