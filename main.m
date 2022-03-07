%% NE
% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 
%                           Game Table:
%  ____________________________________________________________________
% |      W                                    |   W/O                  |
% |___________________________________________|________________________|
% |  W |    u(1,1,1),u(1,1,2)                 |    u(1,2,1),u(1,2,2)   |
% | W/O|    u(2,1,1),u(2,1,2)                 |    u(2,2,1),u(2,2,2)   |
% |____|______________________________________|________________________|
%-------------------------------------------------------------------------

clc
% close all
tic
%% parameters
    N1=2;                                % No.  users
    N2=N1;
    R1=4;                                % BS1 user area center radius
    R2=10;                                % BS2 user area center radius    
    N_i=10;                              % No.  IRS elements
    P_T=(10^(40/10))*1e-3;               % BS power
    [x0,y0,z0]=deal(50,20,0);            % user area center
    [x1,y1,z1]=deal(0,0,0);              % BS1 location
    [x2,y2,z2]=deal(140,0,0);            % BS2 location
    [x_i,y_i,z_i]=deal(50,30,10);        % IRS location
    M_t=1;                               % No. of transmitter antennas
    M_r=1;                               % No. of receiver antennas
    d_IB=100;
    alpha_d=3.76;
    alpha_r=2;
    BW=10e7;
    noise_power=(10^(-114/10));     % -169dbm/Hz
    epsilon=1e-7;
    r=0.23;
    N_iter=200;
    %% users location
    t=2*pi*rand(N1,1);
    radius = 1*sqrt(rand(N1,1));
    x = x0 + radius.*cos(t);
    y = y0 + radius.*sin(t);
    % x(1)=32.52;y(1)=23.48;
    % x(2)=48.45;y(2)=19.55;
    % x(1)=50;y(1)=25;
    % x(2)=40;y(2)=30;
    x(1)=20;y(1)=25;
    x(2)=20;y(2)=25;
    x(3)=50;y(3)=15;
    x(4)=80;y(4)=10;
    x(3)=55;y(3)=10;
    x(4)=55;y(4)=10;
    
    z=1.5*ones(2*N1,1);
    
%%
    utility=zeros(2,2,2,10);
    rate=zeros(2,2,2,10);
    P=zeros(2,10);
    NE_final_rate=zeros(2,10);
    NE_final_utility=zeros(2,10);
    Random_final_rate=zeros(2,10);
    Random_final_utility=zeros(2,10);
    start=8;  
    step=5;    
    m=repmat('.', 1, 1+floor((y_i-start)/step));

    for d=start:step:y_i
        clc
        m(1+floor((d-start)/step))='|';
        disp(m)
        pause(1)
        y(1:2)=d;
        i=1+floor((d-start)/step);
        s=1;                            % s=1 : with IRS
        [utility(:, :, :, i), rate(:, :, :, i)]=game(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter, [R1; R2], r);
        P(:,i)=NE(utility(:, :, :, i));
        NE_final_rate( :, i) = P( 1, i) * P( 2, i) * rate(1 ,1 ,:, i) +...
                                P( 1, i) * ( 1 - P(2, i)) * rate( 1, 2,:, i) + ...
                                (1 - P( 1, i) ) * P( 2, i) * rate( 2, 1, :, i) + ...
                                ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * rate( 2, 2, :, i); 
        NE_final_utility( :, i) = P( 1, i) * P( 2, i) * utility(1 ,1 ,:,i) +...
                               P( 1, i) * ( 1 - P(2, i)) * utility( 1, 2,:, i) + ...
                               (1 - P( 1, i) ) * P( 2, i) * utility( 2, 1, :, i) + ...
                               ( 1 - P( 1, i) ) * (1 - P( 2, i) ) * utility( 2, 2, :, i);

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
    plot((x(1)+x(2))/2,start,'o-','MarkerSize',30,'linewidth',0.1)
    plot((x(1)+x(2))/2,d,'o-','MarkerSize',30,'linewidth',0.1)
    plot((x(1)+x(2))/2,start,'o-','MarkerSize',30,'linewidth',0.1)
    plot((x(1)+x(2))/2,d,'o-','MarkerSize',30,'linewidth',0.1)
%     plot(x(1),y(1),'ro','linewidth',1)
%     plot(x(2),y(2),'bo','linewidth',1)
    plot(x(3),y(3),'ro','linewidth',1)
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
    plot(start:step:y_i,NE_final_rate(1,1:i),'b-*','linewidth',1)
    hold on
    plot(start:step:y_i,Random_final_rate(1,1:i),'b','linewidth',1)
    plot(start:step:y_i,NE_final_rate(2,1:i),'r-*','linewidth',1)
    plot(start:step:y_i,Random_final_rate(2,1:i),'r','linewidth',1)
    grid on
    xlabel('y')
    ylabel('Sum Rate')
    xlim([start d]);
    ylim([0 15]);
    legend('BS1 NE','BS1 RANDOM'...
        ,'BS2 NE','BS2 RANDOM');
    saveas(gcf,'sumrate.jpg')
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
    plot((x(1)+x(2))/2,start,'o-','MarkerSize',30,'linewidth',0.1)
    plot((x(1)+x(2))/2,d,'o-','MarkerSize',30,'linewidth',0.1)
%     plot(x(1),y(1),'ro','linewidth',1)
%     plot(x(2),y(2),'bo','linewidth',1)
    plot(x(3),y(3),'ro','linewidth',1)
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
    plot(start:step:y_i,NE_final_utility(1,1:i),'b-*','linewidth',1)
    hold on
    plot(start:step:y_i,Random_final_utility(1,1:i),'b-o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'linewidth',1)
    plot(start:step:y_i,NE_final_utility(2,1:i),'r-*',...
                'MarkerSize',4,...
                'linewidth',1)
    plot(start:step:y_i,Random_final_utility(2,1:i),'r','linewidth',1)
    grid on
    xlabel('y')
    ylabel('Utility')
    xlim([start d]);
    ylim([0 15]);
    legend('BS1 NE','BS1 RANDOM'...
        ,'BS2 NE','BS2 RANDOM');
    saveas(gcf,'utility.jpg')
    
    