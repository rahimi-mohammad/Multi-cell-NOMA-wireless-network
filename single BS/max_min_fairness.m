clc
close all
%% parameters
N1 = 3 ;                           % No.  users
N_i = 15 ;                         % No.  IRS elements
N_iter = 1 ;                     % No. iteraions
P_T = (10^(40/10)) * 1e-3 ;               % BS power(w)
R = 10 ;                           % Radius
[x0,y0,z0] = deal(30,20,0) ;         % user area center
[x1,y1,z1] = deal(0,0,10) ;       % BS1 location
[x_i,y_i,z_i] = deal(25 * sqrt(2),25 * sqrt(2),10) ;   % IRS location
M_t = 1 ;                          % No. of transmitter antennas
M_r = 1 ;                          % No. of receiver antennas
alpha_d = 3.6 ;
alpha_r = 2 ;
noise_power = (10^(-114/10)) ;     % -169dbm/Hz
interference = 0 ;
epsilon = 0.01 ;
%% users location
t = 2 * pi * rand(N1,1) ;
r  =  R * sqrt(rand(N1,1)) ;
x  =  x0 + r.* cos(t) ;
y  =  y0 + r.* sin(t) ;
z = 1.5 * ones(N1,1) ;
%% 
channel_gain = zeros(N1,1) ;
start = 3 ;
step = 3 ;
final_Lband = zeros(floor(N_i/step) + 1,1) ;
final_rate = zeros(floor(N_i/step) + 1,1) ;
random_final_rate = zeros(floor(N_i/step) + 1,1) ;
tic
for N = start : step : N_i
    m = repmat('.', 1, N_i) ;
    for i = 1 : N
        m(i) = '|' ;m(i + 1) = '|' ;
    end 
    disp(m)
    for k = 1 : N_iter
        %% channel gains
        h_d = (10^0.2) * (randn(N1,M_t,M_r) + 1j * randn(N1,M_t,M_r)) * sqrt(1/2).* ((10^0.2)./(sqrt( (x-x1).^2 + (y-y1).^2 + (z-z1).^2 )).^alpha_d) ;
        g = (10^0.2) * (randn(N,M_t,M_r) + 1j * randn(N,M_t,M_r)) * sqrt(1/2).* ((10^0.2)./(sqrt( (x_i-x1).^2 + (y_i-y1).^2 + (z_i-z1).^2 )).^alpha_r) ;
        h_r = zeros(N,N1) ;
        for i = 1 : N1
            h_r( : ,i) = (10^0.2) * (randn(N,M_t,M_r) + 1j * randn(N,M_t,M_r)) * sqrt(1/2).* ((10^0.2)./(sqrt( (x(i)-x_i).^2 + (y(i)-y_i).^2 + (z(i)-z_i).^2 )).^alpha_r) ;
        end
        %% problem coefficients
        phi = zeros(N,N1) ;
        Q = zeros(N + 1,N + 1,N1) ;
        for i = 1 : N1
            phi( : ,i) = diag(h_r( : ,i)') * g ;
            Q( : , : ,i) = [phi( : ,i) * transpose(conj(phi( : ,i))) conj(h_d(i)) * phi( : ,i) ;h_d(i) * transpose(conj(phi( : ,i))) 0] ;
        end
        %% problem solver
        Q = 10^13 * Q ;
        cvx_begin sdp quiet
            variable X(N + 1,N + 1) complex hermitian   ;
            variable s(1,1)     complex hermitian   ;
            maximize s
            subject to 
            for m = 1 : N1
                10^13 * abs(h_d(m))^2 + real(trace(Q( : , : ,m) * X))>= s
            end
            diag(X) == 1 ;
            X == hermitian_semidefinite(N + 1) ;
        cvx_end
        Z_cvx = GR(X) ;
        Q = 10^(-13) * Q ;
        %% channel gains
        channel_gain = sort(abs(h_d + h_r' * diag(Z_cvx) * g)) ;
        [~, Z_cvx1] = Rate2(h_d, g, h_r, N1, N, P_T, noise_power, N_iter) ;
        Theta = diag(Z_cvx1) ;
        channel_gain1 = sort(abs(h_d + h_r'  *  Theta  *  g)) ;
        [exact_rate, ~] = SISONoma(epsilon, channel_gain1.^2  *  P_T / noise_power) ;
        final_rate(floor(N/step) + 1) = final_rate(floor(N/step) + 1) +  exact_rate ;
        rate_lower_bound = log2(1 + P_T * channel_gain(1)^2/noise_power)/N1 ;
        final_Lband(floor(N/step) + 1) = final_Lband(floor(N/step) + 1) +  rate_lower_bound ;
        %% random theta
        rand_theta = exp(1j * 2 * pi * rand(N,1))  ;
        channel_gain = abs(h_d + h_r' * diag(rand_theta) * g)  ;
        channel_gain = sort(abs(channel_gain))  ;
        rate_lower_bound = log2(1 + P_T * channel_gain(1)^2/noise_power)/N1 ;
        random_final_rate(floor(N/step) + 1) = random_final_rate(floor(N/step) + 1) + rate_lower_bound ;
    
    end
end
    %% without IRS
for k = 1 : N_iter
    h_d = sort(abs(h_d)) ;
    rate_lower_bound = log2(1 + P_T * h_d(1)^2/noise_power)/N1 ;
    final_Lband(1) = final_Lband(1) +  rate_lower_bound ;
    [exact_rate, ~] = SISONoma(epsilon, h_d.^2  *  P_T / noise_power) ;
    final_rate(1) = final_rate(1) + exact_rate ;          
end
final_rate = final_rate/N_iter ;
final_Lband = final_Lband/N_iter ;
random_final_rate = random_final_rate/N_iter ;
toc
%% plot 1 SumRate
figure
tiledlayout(1,2)
nexttile
%  draw BS1
set(gca,'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)  ; 
plot(0,0,'^','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','linewidth',1)
hold on
% draw IRS
plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
% draw BS1 users
plot(x0,y0,'o-','MarkerSize',60  *  R/10,'linewidth',0.1)
% configs
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 150])  ; 
ylim([-5 50])  ; 
legend('BS1','IRS','UE1')  ; 
grid on
nexttile
% draw results
plot([0 step : step : N_i],final_rate, 'g-o','linewidth',1)
hold on
plot([0 step : step : N_i],final_Lband, 'b-o','linewidth',1)
plot([0 step : step : N_i],random_final_rate, 'r-o','linewidth',1)
set(gcf,'color','w');
ylim([0 6])
title('Fair Rate')
xlabel('N')
ylabel('rate(bps/Hz)')
legend('SDR', 'Lower bound','Random')
grid on
saveas(gcf,['rate for ', num2str(N_iter),' Iterations and ', num2str(N1),' users.epsc'])
