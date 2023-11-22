%% NE
% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 
%                           Game Table:
%  ____________________________________________________________________
% |      W/O                                  |     W                  |
% |___________________________________________|________________________|
% | W/O|    u(1,1,1),u(1,1,2)                 |    u(1,2,1),u(1,2,2)   |
% | W  |    u(2,1,1),u(2,1,2)                 |    u(2,2,1),u(2,2,2)   |
% |____|______________________________________|________________________|
%-------------------------------------------------------------------------

clc
close all
tic
%% parameters
    N1=10;                                % No.  users
    N2=5;
    R1=10;                                % BS1 user area center radius
    R2=10;                                % BS2 user area center radius    
    N_i=30;                              % No.  IRS elements
    N_iter=2;
    SV=5;
    P_T=(10^(40/10))*1e-3;               % BS power
%     [x0,y0,z0]=deal(50,20,0);            % user area center
    [x1,y1,z1]=deal(0,0,0);              % BS1 location
    [x2,y2,z2]=deal(100,10,0);            % BS2 location
    [x_i,y_i,z_i]=deal(50,30,10);        % IRS location
    M_t=1;                               % No. of transmitter antennas
    M_r=1;                               % No. of receiver antennas
    d_IB=100;
    alpha_d=3.76;
    alpha_r=2;
    BW=10e7;
    noise_power=(10^(-114/10));     % -169dbm/Hz
    epsilon=1e-7;
%     r=0.23;
    r=0.3;
    %% users location
    x=zeros(N1+N2,1);
    y=zeros(N1+N2,1);
    z=1.5*ones(N1+N2,1);
    x(1:N1,1)=35;
    x(N1+1:N2+N1)=55;
    y(1:N1,1)=25;
    y(N1+1:N2+N1)=10;
%%
    utilityG=zeros(2,2,2,SV);
    rateG=zeros(2,2,2,SV);
    P=zeros(2,SV);
    NE_final_rateG=zeros(2,SV);
    NE_final_utilityG=zeros(2,SV);
    Random_final_rateG=zeros(2,SV);
    Random_final_utilityG=zeros(2,SV);
    start=8;  
    step=5;    
    m=repmat('.', 1, 1+floor((y_i-start)/step));

    for d=start+0*step:step:start+(SV-1)*step
        m(1+floor((d-start)/step))='|';
        disp(m)
        pause(1)
        y(1:N1)=d;
        i=1+floor((d-start)/step);
        s=1;                            % s=1 : with IRS
        [~, rateG(:, :, :, i)]=game(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter, [R1; R2], r);

    %% Random Resource allocation for IRS 
        Random_final_rateG(1, i)=0.5*(rateG(1,2,1, i)+rateG(2,1,1, i));
        Random_final_rateG(2, i)=0.5*(rateG(2,1,2, i)+rateG(1,2,2, i));
    
        Random_final_utilityG(1, i)=0.5*(utilityG(1,2,1, i)+utilityG(2,1,1, i));
        Random_final_utilityG(2, i)=0.5*(utilityG(2,1,2, i)+utilityG(1,2,2, i));

    end
%     calculate utilityG
    for d=start+0*step:step:start+(SV-1)*step
        i=1+floor((d-start)/step);
        utilityG(1, 1, 1,i)=rateG(1, 1, 1,i);
        utilityG(1, 1, 2,i)=rateG(1, 1, 2,i);
        utilityG(1, 2, 1,i)=rateG(1, 2, 1,i);    
        utilityG(2, 1, 2,i)=rateG(2, 1, 2,i);

        utilityG(1, 2, 2,i)=rateG(1, 2, 2,i)-r*N_i;
        utilityG(2, 1, 1,i)=rateG(2, 1, 1,i)-r*N_i;

        utilityG(2, 2, 1,i)=utilityG(2, 1, 1,i)-r*N_i/2;
        utilityG(2, 2, 2,i)=utilityG(1, 2, 2,i)-r*N_i/2;        
        P(:,i)=NE(utilityG(:, :, :, i));
        NE_final_rateG( :, i) = P( 1, i) * P( 2, i) * rateG(1 ,1 ,:, i) +...
                                P( 1, i) * ( 1 - P(2, i)) * rateG( 1, 2,:, i) + ...
                                (1 - P( 1, i) ) * P( 2, i) * rateG( 2, 1, :, i) + ...
                                ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * rateG( 2, 2, :, i); 
        NE_final_utilityG( :, i) = P( 1, i) * P( 2, i) * utilityG(1 ,1 ,:,i) +...
                               P( 1, i) * ( 1 - P(2, i)) * utilityG( 1, 2,:, i) + ...
                               (1 - P( 1, i) ) * P( 2, i) * utilityG( 2, 1, :, i) + ...
                               ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * utilityG( 2, 2, :, i);
    %% Random Resource allocation for IRS 
        Random_final_rate(1, i)=0.5*(rateG(1,2,1, i)+rateG(2,1,1, i));
        Random_final_rate(2, i)=0.5*(rateG(2,1,2, i)+rateG(1,2,2, i));
        Random_final_utility(1, i)=0.5*(utilityG(1,2,1, i)+utilityG(2,1,1, i));
        Random_final_utility(2, i)=0.5*(utilityG(2,1,2, i)+utilityG(1,2,2, i));
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

plot(start:step:(SV-1)*step+start,NE_final_rateG(1,1:i),'b-o','linewidth',1)
hold on
plot(start:step:(SV-1)*step+start,Random_final_rateG(1,1:i),'b','linewidth',1)
plot(start:step:(SV-1)*step+start,NE_final_rateG(2,1:i),'r-o','linewidth',1)
plot(start:step:(SV-1)*step+start,Random_final_rateG(2,1:i),'r','linewidth',1)
grid on
title(['N=',num2str(N_i)])
xlabel('y')
ylabel('Sum rateG')
xlim([start d]);
ylim([-5 10]);
legend('BS1 NE','BS1 RANDOM'...
    ,'BS2 NE','BS2 RANDOM');
saveas(gcf,['rateG for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.epsc'])
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
    % draw results
    plot(start:step:(SV-1)*step+start,NE_final_utilityG(1,1:i),'b-o','linewidth',1)
    hold on
    plot(start:step:(SV-1)*step+start,Random_final_utilityG(1,1:i),'b',...
        'MarkerSize',4,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'linewidth',1)
    plot(start:step:(SV-1)*step+start,NE_final_utilityG(2,1:i),'r-o',...
                'linewidth',1)
    plot(start:step:(SV-1)*step+start,Random_final_utilityG(2,1:i),'r','linewidth',1)
    % configs
    grid on
    title(['N=',num2str(N_i)])
    xlabel('y')
    ylabel('utilityG')
    xlim([start d]);
    ylim([-5 10]);
    legend('BS1 NE','BS1 RANDOM'...
        ,'BS2 NE','BS2 RANDOM');
    saveas(gcf,['utilityG for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.epsc'])
%     