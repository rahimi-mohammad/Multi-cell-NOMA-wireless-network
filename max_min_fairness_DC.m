clc
close all
%% parameters
N1=1;   % No. of users
N=1;    % No. of IRS elements
R=15;   % Radius
[x0,y0,z0]=deal(0.5*25,0,0);  % user area center
[x1,y1,z1]=deal(0,0,0);   % BS1 lozation
[x2,y2,z2]=deal(0.5*100,0,0); % BS2 lozation
[x_i,y_i,z_i]=deal(0.5*50,0.5*30,0.5*10); % IRS lozation
M_t=1;  % No. of transmitter antennas
M_r=1;  % No. of receiver antennas
d_IB=100;
alpha_d=0;
alpha_r=0;
gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
noise_power=0.1;
epsilon=1e-7;
%% users location
t=2*pi*rand(N1,1);
r = R*sqrt(rand(N1,1));
x = x0 + r.*cos(t);
y = y0 + r.*sin(t);
x=0.5*50;
y=0.5*15;
z=zeros(N1,1);
figure
plot(x,y,'ro')
hold on
plot(x_i,y_i,'go')
title('user location')
xlabel('x')
ylabel('y')
channel_gain=zeros(N1,1);

%% channel gains
h_d=zeros(N1,1);    %BS-user link
g=zeros(N,1);       %BS-IRS link
h_r=zeros(N,N1);    %IRS-user link
h_d=(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x-x1).^2+(y-y1).^2+(z-z1).^2 )).^alpha_d);
g=(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^alpha_r);
for i=1:N1
    h_r(:,i)=(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha_r);
end

%% freeze
% N1=1;   % No. of users
% N=2;   % No. of IRS elements
% h_d=zeros(N1,1);
% h_r=zeros(N,N1);
% g=zeros(N,1);
% load('h_d.mat','h_d');
% load('g.mat','g');
% load('h_r.mat','h_r');

%% problem coefficients
phi=zeros(N,N1);
Q=zeros(N+1,N+1,N1);
for i=1:N1
    phi(:,i)=diag(h_r(:,i))*g;
    Q(:,:,i)=[phi(:,i)*phi(:,i)' h_d(i)'*phi(:,i);h_d(i)*phi(:,i)' 0];
end
%% problem solver
cvx_begin sdp
    variable X(N+1,N+1) complex hermitian  ;
    
    maximize real(trace(Q(:,:,1)*X))    %% trace(A4*X) is Real if X and A4 is hemitian matrix
    subject to 
     diag(X) == 1;
     X==hermitian_semidefinite(N+1);
cvx_end
Z_cvx=passive_beamforming(X);
channel_gain=abs(h_d+h_r(:,1)'*diag(Z_cvx)*g);
disp(['CVX channel gain: ',num2str(min(channel_gain))])
disp(['Direct channel gain: ',num2str(min(abs(h_d)))])

%% Analytic Theta for 1-user
analytic_theta=zeros(N,1);
for i=1:N
    analytic_theta(i)=angle(h_d(1))+angle(h_r(i))-angle(g(i));
end
analytic_theta=exp(1j*analytic_theta);
disp(['Analytic channel gain: ',num2str(abs(h_d(1)+h_r'*diag(analytic_theta)*g))])
%% random 
    rand_theta=exp(1j*randi([0,7],N,1));
    disp(['random channel gains:',num2str(abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g))])
    