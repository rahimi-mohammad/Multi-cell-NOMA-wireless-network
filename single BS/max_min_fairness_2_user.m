clc
close all
%% parameters
N1=2;                           % No.  users
N_i=2;                         % No.  IRS elements
N_iter=1;                     % No. iteraions
P_T=(10^(40/10))*1e-3;               % BS power(w)
R=3;                           % Radius
[x0,y0,z0]=deal(10,10,0);         % user area center
[x1,y1,z1]=deal(0,0,0);       % BS1 location
[x2,y2,z2]=deal(0.5*100,0,0);   % BS2 location
[x_i,y_i,z_i]=deal(50,30,10);   % IRS location
M_t=1;                          % No. of transmitter antennas
M_r=1;                          % No. of receiver antennas
% d_IB=100;
alpha_d=3.6;
alpha_r=2;
% gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
% BW=10e7;
noise_power=(10^(-114/10));     % -169dbm/Hz
% epsilon=1e-7;
%% users location
t=2*pi*rand(N1,1);
r = R*sqrt(rand(N1,1));
x = x0 + r.*cos(t);
y = y0 + r.*sin(t);
x(1)=32.52;y(1)=23.48;
x(2)=50;y(2)=15;
% x(1)=50;y(1)=25;
% x(2)=40;y(2)=30;
z=1.5*ones(N1,1);
% y(1)=y(1)-15;
% x(1)=x(1)-20;
%% 
channel_gain=zeros(N1,1);
step=2;
final_rate=zeros(floor(N_i/step)+1,1);
random_final_rate=zeros(floor(N_i/step)+1,1);
tic
for N=N_i:step:N_i
    m=repmat('.', 1, N_i);

    for i=1:N
        m(i)='|';m(i+1)='|';
        
    end
    for k=1:N_iter
        disp(m)
        %% channel gains
        h_d=(10^0.2)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x-x1).^2+(y-y1).^2+(z-z1).^2 )).^alpha_d);
        g=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^alpha_r);
        h_r=zeros(N,N1);
        for i=1:N1
            h_r(:,i)=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha_r);
        end
        %% problem coefficients
        phi=zeros(N,N1);
        Q=zeros(N+1,N+1,N1);
        for i=1:N1
            phi(:,i)=diag(h_r(:,i)')*g;
            Q(:,:,i)=[phi(:,i)*transpose(conj(phi(:,i))) conj(h_d(i))*phi(:,i);h_d(i)*transpose(conj(phi(:,i))) 0];
        end

        %% problem solver
        Q=10^13*Q;
        cvx_begin sdp
            variable X(N+1,N+1) complex hermitian  ;
            maximize min(10^13*abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),10^13*abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
            subject to 
            diag(X) == 1;
            X==hermitian_semidefinite(N+1);
        cvx_end
        Z_cvx=GR(X);
        Q=10^(-13)*Q;
        %% channel gains
        channel_gain(1)=abs(h_d(1)+h_r(:,1)'*diag(Z_cvx)*g);
        channel_gain(2)=abs(h_d(2)+h_r(:,2)'*diag(Z_cvx)*g);
%         channel_gain(1)=sqrt(abs(h_d(1))^2+real(trace(Q(:,:,1)*X)));
%         channel_gain(2)=sqrt(abs(h_d(2))^2+real(trace(Q(:,:,2)*X)));

        rate_lower_bound=(1/N1)*(log2(1+min(channel_gain.^2)*P_T/noise_power))

        channel_gain=sort(abs(channel_gain));
        p = 10^12*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
        r = roots(p);
        alpha=r(0<r&r<1); 
        rate_lower_bound=(1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power))
        final_rate(floor(N/step)+1)=final_rate(floor(N/step)+1)+rate_lower_bound;
        
        
        %% without IRS
%         final_rate(1)=(1/N1)*(log2(1+min(abs(h_d).^2)*P_T/noise_power))+final_rate(1);

        h_d=sort(abs(h_d));
        % h_d=[1;2];
        p = 10^10*[P_T*abs(h_d(1)*h_d(2))^2, noise_power*(h_d'*h_d), -noise_power*abs(h_d(1))^2];
        % p = [P_T*abs(h_d(1)*h_d(2))^2, noise_power*(h_d'*h_d), -noise_power*abs(h_d(2))^2];
        r = roots(p);

        % % if r(1)>0
        alpha=r(0<r&r<1); 
        final_rate(1)=(1/N1)*(log2((abs(h_d(1))^2*P_T+noise_power)/(abs(h_d(1))^2*P_T*alpha+noise_power))+log2((abs(h_d(2))^2*alpha*P_T+noise_power)/noise_power))+final_rate(1);         
        % final_rate(1)=(1/N1)*(log2(1+abs(h_d(1))^2*(alpha)*P_T/noise_power)+log2(1+abs(h_d(2))^2*(1-alpha)*P_T/(noise_power);        
        %% random theta
        rand_theta=exp(1j*2*pi*rand(N,1));
        channel_gain=sort(abs(h_d));
        p = 10^2*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
        r = roots(p);
        alpha=r(0<r&r<1); 
        channel_gain(1)=abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g);
        channel_gain(2)=abs(h_d(2)+h_r(:,2)'*diag(rand_theta)*g);
        rate_lower_bound=(1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power));
% %         rate_lower_bound=(1/N1)*(log2(1+min(channel_gain.^2)*P_T/noise_power));
        random_final_rate(floor(N/step)+1)=random_final_rate(floor(N/step)+1)+rate_lower_bound;
    
        
    end
%     disp('CVX channel gains')
%     
end
final_rate=final_rate/N_iter;
random_final_rate=random_final_rate/N_iter;
toc

%% plot
figure
tiledlayout(1,2)
nexttile

plot(0,0,'^','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
hold on
plot(x(1),y(1),'ro')
hold on
plot(x(2),y(2),'bo')
hold on
plot(x_i,y_i,'gs','MarkerSize',30)
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 100]);
ylim([-5 50]);
legend('BS1','UE1','UE2','IRS');
nexttile
plot([0 step:step:N_i],final_rate, 'r')
hold on
plot([0 step:step:N_i],random_final_rate, 'g')
title('Fair Rate')
xlabel('N')
ylabel('rate(bps/Hz)')
legend('SDR','Random')
