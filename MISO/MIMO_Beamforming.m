    % N1                 - No.  users of BS1,
    % N_i                - No.  IRS elements,
    % P_T                - BS power,
    % x,y,z              - Users location,
    % x1,y1,z1           - BS location,   
    % x_i,y_i,z_i        - IRS location,      
    % alpha_d            - BS to user pathloss
    % alpha_r            - BS to IRS pathloss
    % noise_power        - Noise power
    % N_iter             - No. iteraions

clc
close all
%% parameters
M_t=2;                          % No. of transmitter antennas
M_r=1;                          % No. of receiver antennas
N1=2;                           % No.  users
N_i=3;                         % No.  IRS elements
P_T=(10^(40/10));               % BS power
[x0,y0,z0]=deal(50,20,0);         % user area center
[x1,y1,z1]=deal(100,0,0);       % BS1 location
%     [x2,y2,z2]=deal(0.5*100,0,0);   % BS2 location
[x_i,y_i,z_i]=deal(50,30,10);   % IRS location
alpha_d=3.76;
alpha_r=2;
noise_power=(10^(-114/10));     % -169dbm/Hz
N_iter=1;                      % No. iteraions
step=N_i;
epsilon=1e-7;
gama=100*[1 1 1 1 1 1 1 1];     % minimum SNR
BW=10e7;
d_IB=100;
R=25;                           % Radius
%% users location
x=zeros(N1,1);
y=zeros(N1,1);
z=1.5*ones(N1,1);
x0(1:N1,1)=35;
y0(1:N1,1)=25;
t=2*pi*rand(N1,1);
radius = R(1)*sqrt(rand(N1,1));
x(1:N1) = x0(1:N1) + radius.*cos(t);
y(1:N1) = y0(1:N1) + radius.*sin(t) ;
%% Channel Gain
N=N_i;
[h_d, G, h_r]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                          x_i, y_i, z_i, alpha_d, alpha_r);
                                            
 w=ones(M_t,N1)                                     
 AOA(h_d, G, h_r, w)
                                          %% problem coefficients
phi=zeros(N,N1);
Q=zeros(N+1,N+1,N1);
for i=1:N1
    phi(:,i)=diag(h_r(:,i)')*G;
    Q(:,:,i)=[phi(:,i)*transpose(conj(phi(:,i))) conj(h_d(i))*phi(:,i);h_d(i)*transpose(conj(phi(:,i))) 0];
end

