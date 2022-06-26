% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 
% Usin Potential Game
clc
close all
tic
%% parameters
N1=10;                                % No. first Basestations users
N2=5;                                 % No. first Basestations users
R1=10;                                % BS1 user area center radius
R2=10;                                % BS2 user area center radius    
N_i=6;                               % No.  IRS elements
q_t=1/6;                              % Quantization step-size
N_iter=1;
SV=5;                                 % Scenario Variable
P_T=(10^(40/10))*1e-3;                % BS power
% [x0,y0,z0]=deal(50,20,0);             % user area center
[x1,y1,z1]=deal(0,0,0);               % BS1 location
[x2,y2,z2]=deal(100,10,0);            % BS2 location
[x_i,y_i,z_i]=deal(50,30,10);         % IRS location
M_t=1;                                % No. of transmitter antennas
M_r=1;                                % No. of receiver antennas
d_IB=100;
alpha_d=3.76;
alpha_r=2;
BW=10e7;
noise_power=(10^(-114/10));           % -169dbm/Hz
epsilon=1e-7;
r=0.1;
r=0.3;
%% users location
x=zeros(N1+N2,1);
y=zeros(N1+N2,1);
z=1.5*ones(N1+N2,1);
x(1:N1,1)=35;
x(N1+1:N2+N1)=55;
y(1:N1,1)=25;
y(N1+1:N2+N1)=10;
%% Calcatin Rates and NE
%     utility=zeros(2,2,2,10);
%     rate=zeros(2,2,2,10);
utility=zeros(ceil(1/q_t)+1,ceil(1/q_t)+1,2,SV);
rate=zeros(ceil(1/q_t)+1,ceil(1/q_t)+1,2,SV);
var=[];
save('utility.mat','var')
save('rate.mat','var')
save('loopcounter.mat','var')

P=zeros(2,SV);
NE_final_rate=zeros(2,SV);
NE_final_utility=zeros(2,SV);
Random_final_rate=zeros(2,SV);
Random_final_utility=zeros(2,SV);
step=5;    
start=8;  
m=repmat('.', 1, 1+floor((y_i-start)/step));

for d=start+0*step:step:start+(SV-1)*step
    
    m(1+floor((d-start)/step))='|';
    disp(m)
    pause(1)
    y(1:2)=d;
    
    i=1+floor((d-start)/step);
   [~, rate(:, :, :, i)]=PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                                        x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                                        [R1; R2], r, q_t);
    SaveOut('rate.mat',rate(:,:,:,i)); 
    SaveOut('loopcounter.mat',i)
%% Random Resource allocation for IRS 
%     Random_final_rate(1, i)=0.5*(rate(1,2,1, i)+rate(2,1,1, i));
%     Random_final_rate(2, i)=0.5*(rate(2,1,2, i)+rate(1,2,2, i));
%     Random_final_utility(1, i)=0.5*(utility(1,2,1, i)+utility(2,1,1, i));
%     Random_final_utility(2, i)=0.5*(utility(2,1,2, i)+utility(1,2,2, i));

end

for i=1:SV
    for s1=1:ceil(1/q_t)+1
       for s2=1:ceil(1/q_t)+1
        N_i1=floor((s1-1)*N_i*q_t);    % No. of elements allocated to BS1
        N_i2=floor((s2-1)*N_i*q_t);   % No. of elements allocated to BS2
       utility(s1,s2,1,i)=rate(s1,s2,1,i)-N_i1*r;
       utility(s1,s2,2,i)=rate(s1,s2,2,i)-N_i2*r;
       end            
    end
%     SaveOut('utility.mat',utility(:,:,:,i));
    F=sum(utility(:,:,:,i),3)
    maximum = max(max(F));
    [s1_star,s2_star]=find(F==maximum)
    NE_final_rate( :, i)=rate(s1_star,s2_star,:,i);
    NE_final_utility( :, i)=utility(s1_star,s2_star,:,i);
end
toc
%% plot 1 SumRate
figure
tiledlayout(1,2)
nexttile
% draw BS1
set(gca,'FontName','Times New Roman','FontSize',11,'LineWidth',1.5);
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
plot(x(1),start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
plot(x(1),d,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
% draw BS2 users
plot(x(N1+1),start,'o-','MarkerSize',60*R1/10,'linewidth',0.1)
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
% draw results
set(gca,'FontName','Times New Roman','FontSize',11,'LineWidth',1.5);
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
saveas(gcf,['EPG rate for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.jpg'])
%% plot 2 Utility
figure
tiledlayout(1,2)
nexttile
% draw BS1
set(gca,'FontName','Times New Roman','FontSize',11,'LineWidth',1.5);
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
% draw results
plot(start:step:(SV-1)*step+start,NE_final_utility(1,1:i),'b-*','linewidth',1)
hold on
plot(start:step:(SV-1)*step+start,Random_final_utility(1,1:i),'b',...
    'MarkerSize',4,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_utility(2,1:i),'r-*',...
            'linewidth',1)
plot(start:step:(SV-1)*step+start,Random_final_utility(2,1:i),'r','linewidth',1)
% configs    grid on
title(['N=',num2str(N_i)])
xlabel('y')
ylabel('Utility')
xlim([start max(d, start+step)]);
ylim([0 15]);
legend('BS1 EPG NE','BS1 RANDOM'...
    ,'BS2 EPG NE','BS2 RANDOM');
saveas(gcf,['EPG utility for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.jpg'])
