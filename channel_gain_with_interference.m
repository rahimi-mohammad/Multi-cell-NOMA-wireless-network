% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 
% Average channel gain with interference
%% Before we start
clc
% close all
tic
%% parameters
N1 = 2 ;                                 % No. first Basestations users
N2 = 2 ;                                 % No. first Basestations users
R1 = 10 ;                                % BS1 user area center radius
R2 = 10 ;                                % BS2 user area center radius    
N_i = 24 ;                                % No.  IRS elements
q_t = 1/8 ;                              % Quantization step-size
N_iter = 50 ;
SV = 6 ;                                  % Scenario Variable
P_T = (10^(40/10)) * 1e-3 ;                 % BS power
% [x0 , y0 , z0] = deal(50 , 20 , 0) ;              % user area center
[x1 , y1 , z1] = deal(0 , 0 , 0) ;                % BS1 location
[x2 , y2 , z2] = deal(100 , 10 , 0) ;             % BS2 location
[x_i , y_i , z_i] = deal(50 , 30 , 10) ;          % IRS location
M_t = 1 ;                                 % No. of transmitter antennas
M_r = 1 ;                                 % No. of receiver antennas
d_IB = 100 ; 
alpha_d = 3.76 ; 
alpha_r = 2 ; 
BW = 10e7 ; 
noise_power = (10^(-114/10)) ;            % -169dbm/Hz
r = 0.1 ; 
r = 0.3 ; 
interference = 1 ;                      % Interference between tw BSs.
epsilon = 0.05 ;
N_i1 = floor((3-1) * N_i * q_t) ;     % No. of elements allocated to BS1
N_i2 = floor((4-1) * N_i * q_t) ;     % No. of elements allocated to BS2
%% users location
x = zeros(N1+N2 , 1) ; 
y = zeros(N1+N2 , 1) ; 
z = 1.5 * ones(N1+N2 , 1) ; 
x(1 : N1 , 1) = 35 ; 
x(N1+1 : N2+N1) = 55 ; 
y(1 : N1 , 1) = 25 ; 
y(N1+1 : N2+N1) = 10 ; 
%% Calculating Rates
rate = zeros(ceil(1/q_t)+1 , ceil(1/q_t)+1 , 2 , SV) ; 
var = [] ; 
save('EPG_utility.mat' , 'var')
save('rate.mat' , 'var')
save('loopcounter.mat' , 'var') 
start = 1 ;   
step = 1 ;     
m = repmat('.' ,  1 ,  1 + floor((y_i-start)/step)) ; 
channel_gain11 = zeros( N1 , SV ) ;
channel_gain12 = zeros( N1 , SV) ;
error = zeros( SV , 1) ; % Maximum error
tic
N1 = 2 ;
for N_i2 = floor(N_i*q_t*(start : step : (SV - 1) * step + start))
    m(floor((N_i2 - start)/ step) + 1) = '|' ;
    disp(m)
    i = 1+floor((N_i2/(N_i*q_t)-start)/step) ; 
    rate = zeros(ceil(1/q_t) + 1,ceil(1/q_t) + 1,2) ;
    utility = zeros(ceil(1/q_t) + 1,ceil(1/q_t) + 1,2) ;
    x0 = x ;
    y0 = y ;
    t = 2 * pi * rand(N1,1) ;
    R = [R1 ; R2] ;
    radius  =  R(1) * sqrt(rand(N1,1)) ;
    x(1 : N1) = x0(1 : N1)  +  radius.* cos(t) ;
    y(1 : N1) = y0(1 : N1)  +  radius.* sin(t)  ;

    for j = 1 : N_iter
        if M_t==1
            [h_d1, G1, h_r1] = ChannelGain(M_t, M_r, N1, N_i, x(1 : N1), y(1 : N1), z(1 : N1), x1, y1,...
                                    z1, x_i, y_i, z_i, alpha_d, alpha_r) ;
            [~, Z_cvx1] = Rate2(h_d1, G1(1 : N_i1), h_r1(1 : N_i1, : ), N1, N_i1, P_T, noise_power, N_iter) ;
            %             [~, Z_cvx2] = Rate2(h_d2, G2(1 : N_i2), h_r2(1 : N_i2, : ), N2, N_i2, P_T, noise_power, N_iter) ;
            Z_cvx2 = exp((2 * 1j * pi * rand(N_i2, 1))) ;
            Theta = diag([Z_cvx1 ; interference * Z_cvx2 ;zeros(N_i-N_i1-N_i2,1)]) ;
            channel_gain11( 1 : N1 , i) = channel_gain11(1 : N1 , i) + sort(abs(h_d1 + h_r1' * Theta * G1).^2) ;  % With interference
            interference = 0 ;
            Theta = diag([Z_cvx1 ; interference * Z_cvx2 ;zeros(N_i-N_i1-N_i2,1)]) ;
            channel_gain12(1 : N1 , i) = channel_gain12(1 : N1 , i) + sort(abs(h_d1 + h_r1' * Theta * G1).^2) ; % No interference         
    end
    end 
    
end
channel_gain11 = 1e10 * channel_gain11/N_iter ;
channel_gain12 = 1e10 * channel_gain12/N_iter ;
for N_i2 = floor(N_i*q_t*(start : step : (SV - 1) * step + start))
    i = 1+floor((N_i2/(N_i*q_t)-start)/step) ; 
    error(i, 1) = max(abs((channel_gain12(1 : N1 , i) - channel_gain11(1 : N1 , i) )./...
        channel_gain12(1 : N1 , i) )) ; 
end
figure
plot(2 * ( 1 : SV ),100 * error)
xlabel('Users','FontSize',12,...
       'FontWeight','bold','Color','black')
ylabel('Normalized error(%)','FontSize',12,...
       'FontWeight','bold','Color','black')
set(gcf,'color','w');
grid on
saveas(gcf,['Normalized error for ', num2str(N_iter),' Iterations and at most ', num2str(N1),' users.epsc'])

figure
plot(channel_gain12)
hold on
plot(channel_gain11)

% plot(1e9 * channel_gain12)
    
    
    toc
