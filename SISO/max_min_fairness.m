clc
close all
%% parameters
N1=5;                           % No.  users
N_i=10;                         % No.  IRS elements
N_iter=100;                     % No. iteraions
P_T=(10^(40/10))*1e-3;               % BS power(w)
R=3;                           % Radius
[x0,y0,z0]=deal(10,10,0);         % user area center
[x1,y1,z1]=deal(0,0,10);       % BS1 location
[x_i,y_i,z_i]=deal(25*sqrt(2),25*sqrt(2),10);   % IRS location
M_t=1;                          % No. of transmitter antennas
M_r=1;                          % No. of receiver antennas
alpha_d=3.6;
alpha_r=2;
noise_power=(10^(-114/10));     % -169dbm/Hz
%% users location
t=2*pi*rand(N1,1);
r = R*sqrt(rand(N1,1));
x = x0 + r.*cos(t);
y = y0 + r.*sin(t);
% x(1)=32.52;y(1)=23.48;
x(1)=100;y(1)=5;
x(2)=48.45;y(2)=19.55;
x(2)=0;y(2)=20;
z=1.5*ones(N1,1);
% y(1)=y(1)-15;
% x(1)=x(1)-20;
%% 
channel_gain=zeros(N1,1);
start=2;
step=2;
final_rate=zeros(floor(N_i/step)+1,1);
random_final_rate=zeros(floor(N_i/step)+1,1);
tic
for N=start:step:N_i
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
            variable s(1,1)     complex hermitian  ;
            maximize s
%             min(10^13*abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),10^1e3*abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
            subject to 
            for m=1:N1
                10^13*abs(h_d(m))^2+real(trace(Q(:,:,m)*X))>=s
            end
            diag(X) == 1;
            X==hermitian_semidefinite(N+1);
        cvx_end


        Z_cvx=GR(X);
        Q=10^(-13)*Q;
        %% channel gains
%         for j=1:N1
%             channel_gain(1)=abs(h_d(1)+h_r(:,1)'*diag(Z_cvx)*g);
%             channel_gain(2)=abs(h_d(2)+h_r(:,2)'*diag(Z_cvx)*g);
        channel_gain=sort(abs(h_d+h_r'*diag(Z_cvx)*g));
%         p = 10^12*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
%         r = roots(p);
%         alpha=r(0<r&r<1); 
%         rate_lower_bound=(1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power))
        rate_lower_bound=log2(1+P_T*channel_gain(1)^2/noise_power)/N1;
%         rate_lower_bound=(1/N1)*(log2(1+min(channel_gain.^2)*P_T/noise_power))

        final_rate(floor(N/step)+1)=final_rate(floor(N/step)+1)+rate_lower_bound;
                
        %% without IRS
        h_d=sort(abs(h_d));
%         p = 10^10*[P_T*abs(h_d(1)*h_d(2))^2, noise_power*(h_d'*h_d), -noise_power*abs(h_d(1))^2];
%         r = roots(p);
%         alpha=r(0<r&r<1); 

        final_rate(1)=log2(1+P_T*h_d(1)^2/noise_power)/N1+final_rate(1);         
        %% random theta
        rand_theta=exp(1j*2*pi*rand(N,1));
        channel_gain=abs(h_d+h_r'*diag(rand_theta)*g);
        channel_gain=sort(abs(channel_gain));
%         p = 10^2*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
%         r = roots(p);
%         alpha=r(0<r&r<1); 
%         rate_lower_bound=(1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power));
        rate_lower_bound=log2(1+P_T*channel_gain(1)^2/noise_power)/N1;

        random_final_rate(floor(N/step)+1)=random_final_rate(floor(N/step)+1)+rate_lower_bound;
    
    end
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
    'MarkerFaceColor',[0.5,0.5,0.5],'linewidth',1)
hold on
hold on
plot(x(1),y(1),'ro','linewidth',1)
hold on
plot(x(2),y(2),'bo','linewidth',1)
hold on
plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 100]);
ylim([-5 50]);
grid on
legend('BS1','UE1','UE2','IRS');
nexttile
plot([0 step:step:N_i],final_rate, 'r','linewidth',1)
hold on
plot([0 step:step:N_i],random_final_rate, 'g','linewidth',1)
title('Fair Rate')
xlabel('N')
ylabel('rate(bps/Hz)')
legend('SDR','Random')
grid on
