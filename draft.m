% final comparison
%% plot
figure
tiledlayout(1,2)
nexttile
% draw BS1
plot(0,0,'^','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','linewidth',1)
hold on
% draw BS2
plot(x2,y2,'^','MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','linewidth',1)
% draw IRS
plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
% draw BS1 users
plot((x(1)+x(2))/2,start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(1)+x(2))/2,d,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(1)+x(2))/2,start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(1)+x(2))/2,d,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
%     plot(x(1),y(1),'ro','linewidth',1)
%     plot(x(2),y(2),'bo','linewidth',1)

% draw BS2 users
plot((x(4)+x(5))/2,start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
%     plot(x(4),y(4),'bo','linewidth',1)
% configs
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 150]);
ylim([-5 50]);
legend('BS1','BS2','UE1','UE2','UE1','UE2','IRS');
grid on
nexttile

%% draw results1
plot(start:step:(SV-1)*step+start,NE_final_utilityG(1,1:i),'b-*','linewidth',1)
hold on
plot(start:step:(SV-1)*step+start,Random_final_utility(1,1:i),'b',...
    'MarkerSize',4,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_utilityG(2,1:i),'r-*',...
            'linewidth',1)
plot(start:step:(SV-1)*step+start,Random_final_utility(2,1:i),'r','linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_utility(2,1:i),'r-*',...
            'linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_utility(2,1:i),'b-*',...
            'linewidth',1)
% configs    grid on
title(['N=',num2str(N_i)])
xlabel('y')
ylabel('Utility')
xlim([start max(d, start+step)]);
ylim([0 15]);
legend('BS1 NE','BS1 RANDOM'...
    ,'BS2 NE','BS2 RANDOM');

%% draw results2
% plot(start:step:(SV-1)*step+start,NE_final_utility(1,1:i),'b-*','linewidth',1)
% hold on
% % plot(start:step:(SV-1)*step+start,Random_final_utility(1,1:i),'b',...
% %     'MarkerSize',4,...
% %     'MarkerEdgeColor','b',...
% %     'MarkerFaceColor','b',...
% %     'linewidth',1)
% plot(start:step:(SV-1)*step+start,NE_final_utility(2,1:i),'r-*',...
%             'linewidth',1)
% % plot(start:step:(SV-1)*step+start,Random_final_utility(2,1:i),'r','linewidth',1)
% % configs    grid on
% title(['N=',num2str(N_i)])
% xlabel('y')
% ylabel('Utility')
% xlim([start max(d, start+step)]);
% ylim([0 15]);
% legend('BS1 NE','BS1 RANDOM'...
%     ,'BS2 NE','BS2 RANDOM');
% saveas(gcf,['utility for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.jpg'])




%% plot
figure
tiledlayout(1,2)
nexttile

plot(0,0,'^','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','linewidth',1)
hold on
plot(x2,y2,'^','MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','linewidth',1)
plot((x(1)+x(2))/2,start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(1)+x(2))/2,d,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(1)+x(2))/2,start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(1)+x(2))/2,d,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot((x(4)+x(5))/2,start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
%     plot(x(1),y(1),'ro','linewidth',1)
%     plot(x(2),y(2),'bo','linewidth',1)
%     plot(x(3),y(3),'ro','linewidth',1)
%     plot(x(4),y(4),'bo','linewidth',1)
plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 150]);
ylim([-5 50]);
legend('BS1','BS2','UE1','UE2','UE1','UE2','IRS');
grid on
nexttile
%% draw results1
plot(start:step:(SV-1)*step+start, NE_final_rate(1,1:i),'b-*','linewidth',1)
hold on
plot(start:step:(SV-1)*step+start,Random_final_rate(1,1:i),'b','linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_rate(2,1:i),'r-*','linewidth',1)
plot(start:step:(SV-1)*step+start,Random_final_rate(2,1:i),'r','linewidth',1)
grid on
title(['N=',num2str(N_i)])
xlabel('y')
ylabel('Sum Rate')
xlim([start max(d, start+step)]);
ylim([0 15]);
legend('BS1 NE','BS1 RANDOM'...
    ,'BS2 NE','BS2 RANDOM');
saveas(gcf,['rate for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.jpg'])
%% draw results 2
plot(start:step:(SV-1)*step+start, NE_final_rate(1,1:i),'b-*','linewidth',1)
hold on
plot(start:step:(SV-1)*step+start,Random_final_rate(1,1:i),'b','linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_rate(2,1:i),'r-*','linewidth',1)
plot(start:step:(SV-1)*step+start,Random_final_rate(2,1:i),'r','linewidth',1)
grid on
title(['N=',num2str(N_i)])
xlabel('y')
ylabel('Sum Rate')
xlim([start max(d, start+step)]);
ylim([0 15]);
legend('BS1 NE','BS1 RANDOM'...
    ,'BS2 NE','BS2 RANDOM');
saveas(gcf,['rate for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.jpg'])
