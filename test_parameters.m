clc
close all
%% parameters
N1=2;   % No. of users
N=3;   % No. of IRS elements
R=15;   % Radius
[x0,y0,z0]=deal(0.5*25,0,0);  % user area center
[x1,y1,z1]=deal(0,0,0);   % BS1 lozation
[x2,y2,z2]=deal(0.5*100,0,0); % BS2 lozation
[x_i,y_i,z_i]=deal(0.5*50,0.5*30,0.5*10); % IRS lozation
M_t=1;  % No. of transmitter antennas
M_r=1;  % No. of receiver antennas
d_IB=100;
alpha=2;
alpha2=1;
gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
noise_power=0.1;
epsilon=1e-7;
%% users location
t=2*pi*rand(N1,1);
r = R*sqrt(rand(N1,1));
x = x0 + r.*cos(t);
y = y0 + r.*sin(t);
z=zeros(N1,1);
plot(x,y,'r.')
title('user location')
xlabel('x')
ylabel('y')
%% channel gains
h_d=(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x-x0).^2+(y-y0).^2+(z-z0).^2 )).^alpha);
g=(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x0).^2+(y_i-y0).^2+(z_i-z0).^2 )).^alpha2);
h_r=zeros(N,N1);
for i=1:N1
    h_r(:,i)=(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha2);
end
 %% freeze
% N1=2;   % No. of users
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
    Q(:,:,i)=[phi(:,i)*transpose(conj(phi(:,i))) conj(h_d(i))*phi(:,i);h_d(i)*transpose(conj(phi(:,i))) 0];
end
%% problem solver
cvx_begin sdp
    variable X(N+1,N+1) hermitian complex ;
    
    maximize min(abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
    subject to 
    real(diag(X)) == 1;
    X==hermitian_semidefinite(N+1);
cvx_end
rank(X);
Z_cvx=passive_beamforming(X);
%% channel gains
disp('CVX channel gains')
abs(h_d(1)+h_r(:,1)'*diag(Z_cvx)*g)
abs(h_d(2)+h_r(:,2)'*diag(Z_cvx)*g)
%% random theta
disp('random channel gains')
rand_theta=exp(1j*randi([0,7],N,1));
abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g)
abs(h_d(2)+h_r(:,2)'*diag(rand_theta)*g)
