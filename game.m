clc
close all
%% parameters
    N1=2;                                % No.  users
    N2=N1;
    N_i=20;                              % No.  IRS elements
    P_T=(10^(40/10))*1e-3;               % BS power
    R=10;                                % Radius
    [x0,y0,z0]=deal(50,20,0);            % user area center
    [x1,y1,z1]=deal(0,0,0);              % BS1 location
    [x2,y2,z2]=deal(140,0,0);            % BS2 location
    [x_i,y_i,z_i]=deal(50,30,10);        % IRS location
    M_t=1;                               % No. of transmitter antennas
    M_r=1;                               % No. of receiver antennas
    d_IB=100;
    alpha_d=3.76;
    alpha_r=2;
    gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
    BW=10e7;
    noise_power=(10^(-114/10));     % -169dbm/Hz
    epsilon=1e-7;
    r=0.23;
    %% users location
    t=2*pi*rand(N1,1);
    radius = R*sqrt(rand(N1,1));
    x = x0 + radius.*cos(t);
    y = y0 + radius.*sin(t);
    % x(1)=32.52;y(1)=23.48;
    % x(2)=48.45;y(2)=19.55;
    % x(1)=50;y(1)=25;
    % x(2)=40;y(2)=30;
    x(1)=20;y(1)=25;
    x(2)=25;y(2)=20;
    x(3)=50;y(3)=15;
    x(4)=80;y(4)=10;
    z=1.5*ones(2*N1,1);

%% NE
    
%               Game Table:
%  ____________________________________________________________________
% |      W                                    |   W/O                  |
% |___________________________________________|________________________|
% |  W |    u(1,1,1),u(1,1,2)                 |    u(1,2,1),u(1,2,2)   |
% | W/O|    u(2,1,1),u(2,1,2)                 |    u(2,2,1),u(2,2,2)   |
% |____|______________________________________|________________________|
  
%%
    rates=zeros(2,2,2);
    u=zeros(2,2,2);
    P=zeros(2,10);
    NE_final_rate=zeros(2,10);
    NE_final_utility=zeros(2,10);
    Random_final_rate=zeros(2,10);
    Random_final_utility=zeros(2,10);
    for d=1:10
        
        y(1:2)=[25;20]-d;
        
        s=1;                            % s=1 : with IRS
        rate(1,1,1)=Utility(s,N1,N_i/2,P_T,x(1:2),y(1:2),...
            z(1:2),x1,y1,z1,x_i,y_i,z_i,alpha_d,alpha_r,noise_power);
        rate(1,1,2)=Utility(s,N1,N_i/2,P_T,x(3:4),y(3:4),...
            z(3:4),x2,y2,z2,x_i,y_i,z_i,alpha_d,alpha_r,noise_power);
        rate(1,2,1)=Utility(s,N1,N_i,P_T,x(1:2),y(1:2),...
            z(1:2),x1,y1,z1,x_i,y_i,z_i,alpha_d,alpha_r,noise_power);    
        rate(2,1,2)=Utility(s,N1,N_i,P_T,x(3:4),y(3:4),...
            z(3:4),x2,y2,z2,x_i,y_i,z_i,alpha_d,alpha_r,noise_power);
        s=0;                            % s=0 : without IRS
        rate(1,2,2)=Utility(s,N1,N_i,P_T,x(3:4),y(3:4),...
            z(3:4),x2,y2,z2,x_i,y_i,z_i,alpha_d,alpha_r,noise_power);
        rate(2,1,1)=Utility(s,N1,N_i,P_T,x(1:2),y(1:2),....
            z(1:2),x1,y1,z1,x_i,y_i,z_i,alpha_d,alpha_r,noise_power);
        
        rate(2,2,1)=rate(2,1,1);
        rate(2,2,2)=rate(1,2,2);

        u(1,1,1)=N1*rate(1,1,1)-r*N_i/2;
        u(1,1,2)=N2*rate(1,1,2)-r*N_i/2;
        u(1,2,1)=N1*rate(1,2,1)-r*N_i;    
        u(2,1,2)=N2*rate(2,1,2)-r*N_i;

        u(1,2,2)=N2*rate(1,2,2);
        u(2,1,1)=N1*rate(2,1,1);
        
        u(2,2,1)=u(2,1,1);
        u(2,2,2)=u(1,2,2);
       
        rate(:,:,1)=N1*rate(:,:,1);
        rate(:,:,2)=N2*rate(:,:,2);

        P(:,d)=NE(u);
        NE_final_rate( :, d) = P( 1, d) * P( 2, d) * rate(1 ,1 ,:) +...
                               P( 1, d) * ( 1 - P(2, d)) * rate( 1, 2,:) + ...
                               (1 - P( 1, d) ) * P( 2, d) * rate( 2, 1, :) + ...
                               ( 1 - P( 1, d) ) * (1 - P( 2, d) ) * rate( 2, 2, :);
        NE_final_utility( :, d) = P( 1, d) * P( 2, d) * u(1 ,1 ,:) +...
                                P( 1, d) * ( 1 - P(2, d)) * u( 1, 2,:) + ...
                                (1 - P( 1, d) ) * P( 2, d) * u( 2, 1, :) + ...
                                ( 1 - P( 1, d) ) * (1 - P( 2, d) ) * u( 2, 2, :); 
    %% Random Resource allocation for IRS 
        Random_final_rate(1,d)=0.5*(rate(1,2,1)+rate(2,1,1));
        Random_final_rate(2,d)=0.5*(rate(2,1,2)+rate(1,2,2));
    
        Random_final_utility(1,d)=0.5*(u(1,2,1)+u(2,1,1));
        Random_final_utility(2,d)=0.5*(u(2,1,2)+u(1,2,2));

    end
    %% plot
    tiledlayout(1,2)
    nexttile
    
    plot(0,0,'^','MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b','linewidth',1)
    hold on
    plot(x2,y2,'^','MarkerSize',10,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r','linewidth',1)
    hold on
    plot(x(1),y(1),'ro','linewidth',1)
    hold on
    plot(x(2),y(2),'bo','linewidth',1)
    hold on
    plot(x(3),y(3),'ro','linewidth',1)
    hold on
    plot(x(4),y(4),'bo','linewidth',1)
    plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
    title('User Location')
    xlabel('x')
    ylabel('y')
    xlim([-5 150]);
    ylim([-5 50]);
    legend('BS1','BS2','UE1','UE2','UE1','UE2','IRS');
    grid on
    nexttile
    plot([1:10],NE_final_rate(1,:),'b-*','linewidth',1)
    hold on
    plot([1:10],Random_final_rate(1,:),'b','linewidth',1)
    hold on
    
    plot([1:10],NE_final_rate(2,:),'r-*','linewidth',1)
    hold on

    plot([1:10],Random_final_rate(2,:),'r','linewidth',1)
    hold on
    xlabel('distance')
    ylabel('Sum Rate')
    xlim([0 10]);
    ylim([0 12]);

    legend('BS1 NE-SUM RATE','BS1 RANDOM-SUM RATE'...
        ,'BS2 NE-SUM RATE','BS2 RANDOM-SUM RATE');
    grid on
    saveas(gcf,'sumrate.jpg')