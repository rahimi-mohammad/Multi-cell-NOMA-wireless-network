% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 

clc
close all
tic
%% parameters
N1=10;                                % No. first Basestations users
N2=5;                                 % No. first Basestations users
R1=10;                                % BS1 user area center radius
R2=10;                                % BS2 user area center radius    
N_i=30;                               % No.  IRS elements
q_t=0.2;                              % Quantization step-size
N_iter=100;
SV=2;                                 % Scenario Variable
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
r=0.23;
r=0.6;
%% users location
x=zeros(N1+N2,1);
y=zeros(N1+N2,1);
z=1.5*ones(N1+N2,1);
x(1:N1,1)=35;
x(N1+1:N2+N1)=55;
y(1:N1,1)=25;
y(N1+1:N2+N1)=10;
%%
%     utility=zeros(2,2,2,10);
%     rate=zeros(2,2,2,10);
utility=zeros(ceil(N_i*q_t)+1,ceil(N_i*q_t)+1,2,SV);
rate=zeros(ceil(N_i*q_t)+1,ceil(N_i*q_t)+1,2,SV);
% var=[];
% save('utility.mat','var')
% save('rate.mat','var')
P=zeros(2,SV);
NE_final_rate=zeros(2,SV);
NE_final_utility=zeros(2,SV);
Random_final_rate=zeros(2,SV);
Random_final_utility=zeros(2,SV);
step=5;    
start=8+3*step;  
m=repmat('.', 1, 1+floor((y_i-start)/step));

for d=start+0*step:step:start+(SV-1)*step
    clc
    m(1+floor((d-start)/step))='|';
    disp(m)
    pause(1)
    y(1:2)=d;
    
    i=1+floor((d-start)/step);
%         [utility(:, :, :, i), rate(:, :, :, i)]=game(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
%                 x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter, [R1; R2], r);
%         P(:,i)=NE(utility(:, :, :, i));
%         NE_final_rate( :, i) = P( 1, i) * P( 2, i) * rate(1 ,1 ,:, i) +...
%                                 P( 1, i) * ( 1 - P(2, i)) * rate( 1, 2,:, i) + ...
%                                 (1 - P( 1, i) ) * P( 2, i) * rate( 2, 1, :, i) + ...
%                                 ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * rate( 2, 2, :, i); 
%         NE_final_utility( :, i) = P( 1, i) * P( 2, i) * utility(1 ,1 ,:,i) +...
%                                P( 1, i) * ( 1 - P(2, i)) * utility( 1, 2,:, i) + ...
%                                (1 - P( 1, i) ) * P( 2, i) * utility( 2, 1, :, i) + ...
%                                ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * utility( 2, 2, :, i);


%         P(:,i)=NE(utility(:, :, :, i));
%         NE_final_rate( :, i) = P( 1, i) * P( 2, i) * rate(1 ,1 ,:, i) +...
%                                 P( 1, i) * ( 1 - P(2, i)) * rate( 1, 2,:, i) + ...
%                                 (1 - P( 1, i) ) * P( 2, i) * rate( 2, 1, :, i) + ...
%                                 ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * rate( 2, 2, :, i); 
%         NE_final_utility( :, i) = P( 1, i) * P( 2, i) * utility(1 ,1 ,:,i) +...
%                                P( 1, i) * ( 1 - P(2, i)) * utility( 1, 2,:, i) + ...
%                                (1 - P( 1, i) ) * P( 2, i) * utility( 2, 1, :, i) + ...
%                                ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * utility( 2, 2, :, i);

   [utility(:, :, :, i), rate(:, :, :, i)]=PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                                        x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                                        [R1; R2], r, q_t);
    SaveOut('utility.mat',utility(:,:,:,i));
    SaveOut('rate.mat',rate(:,:,:,i));                                    
    F=sum(utility,3);
    maximum = max(max(F));
    [s1_star,s2_star]=find(F==maximum);
%     utility(x,y,:,:)
%% Random Resource allocation for IRS 
    Random_final_rate(1, i)=0.5*(rate(1,2,1, i)+rate(2,1,1, i));
    Random_final_rate(2, i)=0.5*(rate(2,1,2, i)+rate(1,2,2, i));
    Random_final_utility(1, i)=0.5*(utility(1,2,1, i)+utility(2,1,1, i));
    Random_final_utility(2, i)=0.5*(utility(2,1,2, i)+utility(1,2,2, i));

end

toc
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
legend('BS1 NE','BS1 RANDOM'...
    ,'BS2 NE','BS2 RANDOM');
saveas(gcf,['utility for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.jpg'])
