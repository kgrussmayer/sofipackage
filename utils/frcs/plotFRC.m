%% plot FRC
clear frc;
clear legendinfo;

io=2;
for ii = 1:size(results,1)
    for jj=1:7
   frc(ii,jj)=results(ii,io).resMid(jj);
   legendinfo{jj} = ['Sector: ',num2str(jj)];
    end
    
end

figure, 
plot(1:length(frc),frc,'x--');
xlabel('Frame of SR sequence');
ylabel('sFRC');
legend(legendinfo);

% ylim([100 300]);
% xlim([0.5 4.5]);
set(gca,'XTick',[1:7])
% xframes = [1, 2, 3, 4,];
% set(gca,'XTickLabel',xframes);
text(1:length(frc),frc,num2str(round(10*frc)./10),'FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','bottom')