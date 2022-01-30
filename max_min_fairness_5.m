clc
close all
%% parameters
N1=2;   % No.  users
N_i=30;   % No.  IRS elements
N_iter=1;  % No. iteraions
P_T=1;    % BS power
R=25;   % Radius
[x0,y0,z0]=deal(0,0,0);  % user area center
[x1,y1,z1]=deal(100,0,0);   % BS1 lozation
[x2,y2,z2]=deal(0.5*100,0,0); % BS2 lozation
[x_i,y_i,z_i]=deal(100,50,0); % IRS lozation
M_t=1;  % No. of transmitter antennas
M_r=1;  % No. of receiver antennas
d_IB=100;
alpha=2;
alpha2=2;
gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
BW=10e7;
noise_power=(10^(-114/10));            % -169dbm/Hz
epsilon=1e-7;
%% users location
t=2*pi*rand(N1,1);
r = R*sqrt(rand(N1,1));
x = x0 + r.*cos(t);
y = y0 + r.*sin(t);
z=zeros(N1,1);
figure
plot(x(1),y(1),'ro')
hold on
plot(x(2),y(2),'bo')
hold on
plot(x_i,y_i,'go')
title('user location')
xlabel('x')
ylabel('y')
channel_gain=zeros(N1,1);
final_rate=zeros(N_i+1,1);
for N=1:1:N_i
    for k=1:N_iter
        %% channel gains
        h_d=(10^0.2)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x-x0).^2+(y-y0).^2+(z-z0).^2 )).^alpha);
        g=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x0).^2+(y_i-y0).^2+(z_i-z0).^2 )).^alpha2);
        h_r=zeros(N,N1);
        for i=1:N1
            h_r(:,i)=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha2);
        end
        %% problem coefficients
        phi=zeros(N,N1);
        Q=zeros(N+1,N+1,N1);
        for i=1:N1
            phi(:,i)=diag(h_r(:,i))*g;
            Q(:,:,i)=[phi(:,i)*transpose(conj(phi(:,i))) conj(h_d(i))*phi(:,i);h_d(i)*transpose(conj(phi(:,i))) 0];
        end
        %% problem solver
        cvx_begin sdp
            variable X(N+1,N+1) complex hermitian  ;
            maximize min(abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
            subject to 
            diag(X) == 1;
            X==hermitian_semidefinite(N+1);
        cvx_end
        rank(X);
        Z_cvx=passive_beamforming(X);
        %% channel gains
        channel_gain(1)=abs(h_d(1)+h_r(:,1)'*diag(Z_cvx)*g);
        channel_gain(2)=abs(h_d(2)+h_r(:,2)'*diag(Z_cvx)*g);
        rate_lower_bound=(1/N1)*(log2(1+min(channel_gain.^2)*P_T/noise_power));
        final_rate(N+1)=final_rate(N+1)+rate_lower_bound;
    end
    disp('CVX channel gains')

%     %% random theta
%     disp('random channel gains')
%     rand_theta=exp(1j*randi([0,7],N,1));
%     abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g);
%     abs(h_d(2)+h_r(:,2)'*diag(rand_theta)*g);
%     disp('difference:')
%     channel_gain-[abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g);abs(h_d(2)+h_r(:,2)'*diag(rand_theta)*g)];
end
final_rate=final_rate/N_iter;
%% without IRS
final_rate(1)=(1/N1)*(log2(1+min(abs(h_d).^2)*P_T/noise_power));
%% plot
figure
plot([0 1:1:N_i],final_rate)
xlabel('N')
ylabel('rate(bps/Hz)')
