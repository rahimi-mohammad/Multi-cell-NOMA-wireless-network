%% NE
% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% main.m 
% This script generates data for train and test of a ML model for solving
% our Optimization problem
tic
%% parameters
N_data=3;
N1=5;                           % No.  users
N_i=10;                         % No.  IRS elements
% N_iter=100;                     % No. iteraions
P_T=(10^(40/10))*1e-3;               % BS power(w)
R=10;                           % Radius
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
z=1.5*ones(N1,1);
filename = 'testdata.xlsx';
xlRange='A1';
for j=1:N_data
    [h_d, G, h_r]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                    x_i, y_i, z_i, alpha_d, alpha_r);
     Data=[transpose(h_d); G; reshape(h_r,[],1)]; 
    final_rate=0;
    N=N_i;
    
    if N_i==0
        % without IRS
                if M_t==1
                    %% without IRS
                    h_d=sort(abs(h_d));
                final_rate=log2(1+P_T*h_d(1)^2/noise_power)/N1+final_rate; 
%                     final_rate=(1/N1)*(log2((abs(h_d(1))^2*P_T+noise_power)/(abs(h_d(1))^2*P_T*alpha+noise_power))+log2((abs(h_d(2))^2*alpha*P_T+noise_power)/noise_power))+final_rate;         
                elseif M_t>1
                    % I will write this later
                    [h_d, ~, ~]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                            x_i, y_i, z_i, alpha_d, alpha_r);      
                end 
    else
        % with IRS
            
            if M_t==1
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
                        subject to 
                        for m=1:N1
                            10^13*abs(h_d(m))^2+real(trace(Q(:,:,m)*X))>=s
                        end
                        diag(X) == 1;
                        X==hermitian_semidefinite(N+1);
                    cvx_end
                    Z_cvx=GR(X);
    %                 Q=10^(-13)*Q;
                    %% channel gains
                    channel_gain=sort(abs(h_d+h_r'*diag(Z_cvx)*g));
                    rate_lower_bound=log2(1+P_T*channel_gain(1)^2/noise_power)/N1;
                    final_rate=final_rate+rate_lower_bound;
            elseif M_t>1
                % I will write this later
                    %% problem coefficients
                    phi=zeros(N,N1);
                    Q=zeros(N+1,N+1,N1);
                    for i=1:N1
                        phi(:,i)=diag(h_r(:,i)')*G;
                        Q(:,:,i)=[phi(:,i)*transpose(conj(phi(:,i))) conj(h_d(i))*phi(:,i);h_d(i)*transpose(conj(phi(:,i))) 0];
                    end
            end
    end
    
    rate=final_rate;
    xlRange=['A', num2str((j-1)*(N1+N_i+N1*N_i)+1)]
    A = {...
        num2cell(real(Data)), num2cell(imag(Data)),...
        num2cell([angle(Z_cvx); zeros(N1+N_i+N1*N_i-N_i, 1) ]),...
        num2cell([rate; zeros(N1+N_i+N1*N_i-1,1)])...
    }
    xlswrite(filename,A,sheet,xlRange)
end









