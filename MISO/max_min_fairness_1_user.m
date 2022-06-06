clc
close all
%% parameters
    N1=1;                           % No.  users
    N_i=30;                         % No.  IRS elements
    N_iter=1;                      % No. iteraions
    P_T=(10^(40/10));               % BS power
    R=25;                           % Radius
    [x0,y0,z0]=deal(50,20,0);         % user area center
    [x1,y1,z1]=deal(100,0,0);       % BS1 location
    [x2,y2,z2]=deal(0.5*100,0,0);   % BS2 location
    [x_i,y_i,z_i]=deal(50,30,10);   % IRS location
    M_t=1;                          % No. of transmitter antennas
    M_r=1;                          % No. of receiver antennas
    d_IB=100;
    alpha_d=3.76;
    alpha_r=2;
    gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
    BW=10e7;
    noise_power=(10^(-114/10));     % -169dbm/Hz
    epsilon=1e-7;
    step=N_i;
    final_rate=zeros(floor(N_i/step)+1+1,1);

%% users location
    t=2*pi*rand(N1,1);
    r = R*sqrt(rand(N1,1));
    x = x0 + r.*cos(t);
    y = y0 + r.*sin(t);
    x(1)=40;y(1)=25;
    z=zeros(N1,1);
    figure
    plot(0,0,'^','MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0.5,0.5,0.5])
    hold on
    plot(x,y,'ro')
    hold on
    plot(x_i,y_i,'gs','MarkerSize',30)
    title('user location')
    xlabel('x')
    ylabel('y')
    xlim([-5 100]);
    ylim([-5 50]);
    legend('BS','UE1','IRS');
    channel_gain=zeros(N1,1);
%%
channel_gain=zeros(N1,1);
step=1;
final_rate=zeros(floor(N_i/step)+1,1);
analytic_final_rate=zeros(floor(N_i/step)+1,1);

for N=1:step:N_i    
%% channel gains
    h_d=zeros(N1,1);    %BS-user link
    g=zeros(N,1);       %BS-IRS link
    h_r=zeros(N,N1);    %IRS-user link
    h_d=(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x-x1).^2+(y-y1).^2+(z-z1).^2 )).^alpha_d);
    g=(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^alpha_r);
    for i=1:N1
        h_r(:,i)=(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha_r);
    end

%% problem coefficients
    phi=zeros(N,N1);
    Q=zeros(N+1,N+1,N1);
    for i=1:N1
        phi(:,i)=diag(h_r(:,i)')*g;
        Q(:,:,i)=[phi(:,i)*phi(:,i)' h_d(i)'*phi(:,i);h_d(i)*phi(:,i)' 0];
    end
    Q=10^14*Q;
    
%% problem solver
    cvx_begin sdp
        variable X(N+1,N+1) complex hermitian  ;

        maximize real(trace(Q(:,:,1)*X))    %% trace(A4*X) is Real if X and A4 is hemitian matrix
        subject to 
         diag(X) == 1;
         X==hermitian_semidefinite(N+1);
    cvx_end
    Z_cvx=GR(X);
    channel_gain=abs(h_d+h_r(:,1)'*diag(Z_cvx)*g)^2;
    Q=10^(-14)*Q;
    final_rate(floor(N/step)+1)=log2((channel_gain*P_T+noise_power)/noise_power);

%% Analytic Theta for 1-user
    analytic_theta=zeros(N,1);
    for i=1:N
        analytic_theta(i)=angle(h_d(1))+angle(h_r(i))-angle(g(i));
    end
    analytic_theta=exp(1j*analytic_theta);
    analytic_final_rate(floor(N/step)+1)=log2((abs(h_d(1)+h_r(:,1)'*diag(analytic_theta)*g)^2*P_T+noise_power)/noise_power);
end
    
    %% Without RIS
    final_rate(1)=log2((abs(h_d(1))^2*P_T+noise_power)/noise_power);
%% random 
%     rand_theta=exp(1j*2*pi*rand(N,1));
%     disp(['random channel gain:',num2str(abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g)^2)])
%     final_rate(2)=log2((abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g)^2*P_T+noise_power)/noise_power);
%% comparisopn
%     analytic_Z=[analytic_theta;1]*[analytic_theta;1 ]';
%     disp(['|h_d|^2: ',num2str(abs(h_d)^2)])
%     disp(['|h_d+h_r*Z_cvx*g|^2: ',num2str(channel_gain)])
%     disp(['rate: ',num2str(final_rate')]);
%     disp(['trace(Q*X)+|h_d|^2=  ',num2str(real(trace(Q(:,:,1)*X))+abs(h_d(1))^2)])
%     disp(['|h_d+h_r*analytic_theta*g|^2=  ',num2str(abs(h_d(1)+h_r'*diag(analytic_theta)*g)^2)])
%     disp(['',num2str(analytic_Z)]);
%% plot
figure
tiledlayout(1,2)
nexttile

plot(0,0,'^','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
hold on
plot(x(1),y(1),'ro','linewidth',1)
hold on
plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 100]);
ylim([-5 50]);
legend('BS1','UE1','IRS');
nexttile
plot(step:step:N_i,final_rate(2:end), 'r','linewidth',1)
hold on
plot(step:step:N_i,analytic_final_rate(2:end), 'g-*')
title('Rate')
xlabel('N')
ylabel('rate(bps/Hz)')
% xlim([1 40]);
% ylim([-5 50]);
legend('SDR','Analytical')
saveas(gcf,'SISO-one-user.jpg')
