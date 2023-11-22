%% Calculating the exact rate of users with max-min fairness objective
clc
close all
%% parameters
N1 = 4;                           % No.  users
N2 = 1;
N_i = 8;                         % No.  IRS elements
N_iter = 1;                     % No. iteraions
P_T = (10^(30/10))*1e-3;               % BS power(w)
R = 10;                           % Radius
[x0,y0,z0] = deal(10,10,0);         % user area center
[x1,y1,z1] = deal(40,0,0);       % BS1 location
[x_i,y_i,z_i] = deal(100,50,0) ;          % IRS location
M_t = 1;                          % No. of transmitter antennas
M_r = 1;                          % No. of receiver antennas
% d_IB = 100 ;
alpha_d = 3.6;
alpha_r = 2;
% alpha_r2 = 2.5;
% BW = 10e7;
% noise_power = 0.01 ;
noise_power = (10^(-114/10));     % -169dbmHz
% epsilon = 1e-7;
interference = 0 ;                      % Interference between tw BSs.
%% users location
x  =  zeros(N1+N2,1) ; 
y  =  zeros(N1+N2,1) ; 
z  =  1.5 * ones(N1+N2,1) ; 
x(1:N1,1)  =  0 ; 
x(N1+1:N2+N1)  =  55 ; 
y(1:N1,1)  =  0 ; 
y(N1+1:N2+N1)  =  10 ; 
x0 = x ;
y0 = y ;
t = 2 * pi * rand(N1,1) ;
radius  =  R(1) * sqrt(rand(N1,1)) ;
x(1 : N1) = x0(1 : N1)  +  radius.* cos(t) ;
y(1 : N1) = y0(1 : N1)  +  radius.* sin(t)  ;
[h_d1, G1, h_r1] = ChannelGain(M_t, M_r, N1, N_i, x(1 : N1), y(1 : N1), z(1 : N1), x1, y1,...
                                        z1, x_i, y_i, z_i, alpha_d, alpha_r) ;
[~, Z_cvx1] = Rate2(h_d1, G1(1 : N_i), h_r1(1 : N_i, : ), N1, N_i, P_T, noise_power, N_iter) ;
Theta = diag(Z_cvx1 ) ;
channel_gain = sort(abs(h_d1 + h_r1' * Theta * G1), "ascend") ;
% channel_gain1 = sqrt([1 ; 2; 3] / 3)  ;
channel_gain1 = channel_gain.^2 * P_T / noise_power ; % normalized channel gain
%% calculating optimal rate
stepsize = 0.01 ;
q_min = 1 ;
q_max = 50 ;
b_star = zeros( N1+1, 0) ;
epsilon = 0.05 ;
while q_max - q_min >= epsilon
    q = (q_min + q_max)/2 ;
    cvx_begin quiet
        dual variable y
        variable b(N1 + 1)   ;
        minimize -q 
        subject to 
        b(1) == 1 ;
        b(N1 + 1) == 0 ;
        for m = 1 : N1 
            b(m + 1) - b(m) <= 0
            -(1 + b(m) * channel_gain1(m)) + q + b(m + 1) * channel_gain1(m) * q <= 0
        end
    cvx_end
    if cvx_status( 1 : 6 ) == 'Infeas'
        q_max = q ;
    else
        q_min = q ;
        b_star = b ;
    end
end
disp(['Q = ', num2str(log2(q))])
for m = 1 : N1 
    disp(['R_', num2str(m),' = ',...
        num2str(log2((1 + channel_gain1(m) * b_star(m))/(1 + channel_gain1(m) * b_star(m + 1)))) ...
        ])
end
rate_lower_bound1 = log2(1 + P_T * channel_gain(1)^2 / noise_power) / N1 ;
disp(['rate_lower_bound1 = ', num2str(rate_lower_bound1)])
b_star(2)-((1/q-1)/channel_gain1(1)+1/q)
Q = q ;
q = channel_gain1(1)^(1/(N1 - 1)) ;
b(1) = 1 ;
b(N1 + 1) = 0 ;
for k = 1 : N1 - 1
    b(k + 1) = (1/q - 1)/ channel_gain1(k) + b(k)/q ;
end
b > 0
b_star = b;
for m = 1 : N1 
    disp(['R_', num2str(m),' = ',...
        num2str(log2((1 + channel_gain1(m) * b_star(m))/(1 + channel_gain1(m) * b_star(m + 1)))) ...
        ])
end
%% 
% channel_gain = zeros(N1,1);
% step = 2;
% final_rate = zeros(floor(N_i/step)+1,1);
% random_final_rate = zeros(floor(N_i/step)+1,1);
% tic
% for N = N_i:step:N_i
%     m = repmat('.', 1, N_i);
% 
%     for i = 1:N
%         m(i) = '|';m(i+1) = '|';
%         
%     end
%     for k = 1:N_iter
%         disp(m)
%         %% channel gains
%         h_d = (10^-1.5)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r)).*...
%             (sqrt(1/2)./(sqrt( (x-x1).^2+(y-y1).^2+(z-z1).^2 )).^(alpha_d/2));
%         g = (10^-1.5)*( sqrt(k1/(1+k1)) +sqrt(1/(1+k1))*(randn(N,M_t,M_r)+ 1j* randn(N,M_t,M_r)) )*...
%             sqrt(1/2).*(1./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^(alpha_r1/2));
%         h_r = zeros(N,N1);
%         for i = 1:N1
%             h_r(:,i) = (10^-1.5)*(sqrt(k2/(1+k2)) +sqrt(1/(1+k2))*...
%                     (randn(N,M_t,M_r)+1j*randn(N,M_t,M_r)))*...
%                     sqrt(1/2).*(1./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^(alpha_r2/2));
%         end
%         %% problem coefficients
%         phi = zeros(N,N1);
%         Q = zeros(N+1,N+1,N1);
%         for i = 1:N1
%             phi(:,i) = diag(h_r(:,i)')*g;
%             Q(:,:,i) = [phi(:,i)*transpose(conj(phi(:,i))) conj(h_d(i))*phi(:,i);h_d(i)*transpose(conj(phi(:,i))) 0];
%         end
% 
%         %% problem solver
%         Q = 10^13*Q;
%         cvx_begin sdp
%             variable X(N+1,N+1) complex hermitian  ;
%             maximize min(10^13*abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),10^13*abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
%             subject to 
%             diag(X)  ==  1;
%             X == hermitian_semidefinite(N+1);
%         cvx_end
%         Z_cvx = GR(X);
%         Q = 10^(-13)*Q;
%         %% channel gains
%         channel_gain(1) = abs(h_d(1)+h_r(:,1)'*diag(Z_cvx)*g);
%         channel_gain(2) = abs(h_d(2)+h_r(:,2)'*diag(Z_cvx)*g);
% %         channel_gain(1) = sqrt(abs(h_d(1))^2+real(trace(Q(:,:,1)*X)));
% %         channel_gain(2) = sqrt(abs(h_d(2))^2+real(trace(Q(:,:,2)*X)));
% 
% %         rate_lower_bound = (1/N1)*(log2(1+min(channel_gain.^2)*P_T/noise_power));
% 
%         channel_gain = sort(abs(channel_gain));
%         p  =  10^12*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
%         r  =  roots(p);
%         alpha = r(0<r&r<1); 
%         rate_lower_bound = (1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power));
%         final_rate(floor(N/step)+1) = final_rate(floor(N/step)+1)+rate_lower_bound;
%                 
%         %% without IRS
%         h_d = sort(abs(h_d));
%         p  =  10^10*[P_T*abs(h_d(1)*h_d(2))^2, noise_power*(h_d'*h_d), -noise_power*abs(h_d(1))^2];
%         r  =  roots(p);
%         alpha = r(0<r&r<1); 
%         final_rate(1) = (1/N1)*(log2((abs(h_d(1))^2*P_T+noise_power)/(abs(h_d(1))^2*P_T*alpha+noise_power))+log2((abs(h_d(2))^2*alpha*P_T+noise_power)/noise_power))+final_rate(1);         
%         %% random theta
%         rand_theta = exp(1j*2*pi*rand(N,1));
%         channel_gain(1) = abs(h_d(1)+h_r(:,1)'*diag(rand_theta)*g);
%         channel_gain(2) = abs(h_d(2)+h_r(:,2)'*diag(rand_theta)*g);
%         channel_gain = sort(abs(channel_gain));
%         p  =  10^2*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
%         r  =  roots(p);
%         alpha = r(0<r&r<1); 
%         rate_lower_bound = (1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power));
%         random_final_rate(floor(N/step)+1) = random_final_rate(floor(N/step)+1)+rate_lower_bound;
%     
%         
%     end
%   
% end
% final_rate = final_rate/N_iter;
% random_final_rate = random_final_rate/N_iter;
% toc
% 
% %% plot
% figure
% tiledlayout(1,2)
% nexttile
% plot(0,0,'^','MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5],'linewidth',1)
% hold on
% hold on
% plot(x(1),y(1),'ro','linewidth',1)
% hold on
% plot(x(2),y(2),'bo','linewidth',1)
% hold on
% plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
% title('User Location')
% xlabel('x')
% ylabel('y')
% xlim([-5 100]);
% ylim([-5 50]);
% grid on
% legend('BS1','UE1','UE2','IRS');
% nexttile
% plot([0 step:step:N_i],final_rate, 'r','linewidth',1)
% hold on
% plot([0 step:step:N_i],random_final_rate, 'g','linewidth',1)
% title('Fair Rate')
% xlabel('N')
% ylabel('rate(bps/Hz)')
% legend('SDR','Random')
% grid on
