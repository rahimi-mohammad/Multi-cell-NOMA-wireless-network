% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 
% Using Potential Game
%% Before we start
clc
close all
tic
%% parameters
N1 = 5 ;                                 % No. first Basestations users
N2 = 1 ;                                 % No. first Basestations users
R1 = 10 ;                                % BS1 user area center radius
R2 = 10 ;                                % BS2 user area center radius    
N = 6 ;                                % No.  IRS elements
q_t = 1/8 ;                              % Quantization step-size
N_iter = 10 ;
SV = 1 ;                                  % Scenario Variable
P_T = (10^(40/10)) * 1e-3 ;                 % BS power
% [x0 , y0 , z0] = deal(50 , 20 , 0) ;              % user area center
[x1 , y1 , z1] = deal(0 , 0 , 0) ;                % BS1 location
[x2 , y2 , z2] = deal(100 , 10 , 0) ;             % BS2 location
[x_i , y_i , z_i] = deal(50 , 30 , 10) ;          % IRS location
M_t = 4 ;                                 % No. of transmitter antennas
M_r = 1 ;                                 % No. of receiver antennas
d_IB = 100 ; 
alpha_d = 3.76 ; 
alpha_r = 2 ; 
% alpha_d = 0 ; 
% alpha_r = 0 ; 
BW = 10e7 ; 
noise_power = (10^(-114/10)) ;            % -169dbm/Hz
interference = 0 ;                      % Interference between two BSs.
epsilon = 0.05 ;
%% users center location
x0 = zeros(N1+N2 , 1) ; 
y0 = zeros(N1+N2 , 1) ; 
z = 1.5 * ones(N1+N2 , 1) ; 
x0(1 : N1 , 1) = 35 ; 
x0(N1+1 : N2+N1) = 55 ; 
y0(1 : N1 , 1) = 25 ; 
y0(N1+1 : N2+N1) = 10 ;
R = [10 ; 10] ;
%%
start = 1 ;
step = 1 ;
final_rate = zeros(SV , 1) ;
tic
for d = start + 0 * step : step : start + (SV-1) * step
%     N_iter = d * 100 ;
    for i = 1 : N_iter
        m = repmat('|' , 1 , i) ;
        disp(m)
    t = 2 * pi * rand(N1,1) ;
    radius  =  R(1) * sqrt(rand(N1,1)) ;
    x(1 : N1 , 1) = x0(1 : N1)  +  radius.* cos(t) ;
    y(1 : N1 , 1) = y0(1 : N1)  +  radius.* sin(t)  ;
    [h_d1, G1, h_r1] = ChannelGain(M_t, M_r, N1, N, x(1 : N1 , 1), y(1 : N1 , 1), z(1 : N1), x1, y1,...
                                        z1, x_i, y_i, z_i, alpha_d, alpha_r) ;
    [W , rate , Z_cvx] = Rate3(h_d1, G1, h_r1 , N1, N, P_T, noise_power , epsilon) ;
    final_rate(d) = final_rate(d) + rate ;
    end
end
final_rate = final_rate/N_iter ;
toc

