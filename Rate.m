% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method calculates aceivable rate of the BS's users for each
% strategy
% Inputs:
%     s                  - Scenario: s=0 without IRS and s=1 with RIS,
    % N1                 - No.  users of BS1,
    % N_i                - No.  IRS elements(It can be a vector),
    % P_T                - BS power,
    % x,y,z              - Users location,
    % x1,y1,z1           - BS location,   
    % x_i,y_i,z_i        - IRS location,      
    % alpha_d            - BS to user pathloss
    % alpha_r            - BS to IRS pathloss
    % noise_power        - Noise power
    % N_iter             - No. iteraions
% Outputs:
    % U                  - Acheivable Rate of the BS(the lower bound).
% ------------------------------------------------------------------------
function rate=Rate(s, M_t, M_r, N1, N_i, P_T, x, y, z, x1, y1, z1,...
                    x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter)
    %% parameters
    channel_gain=zeros(N1,1);
    final_rate=0;
    N=N_i;
    switch s
        case{0} % without IRS
            for k=1:N_iter
                %% channel gains
                [h_d, ~, ~]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                        x_i, y_i, z_i, alpha_d, alpha_r);               
                if M_t==1
                    %% without IRS
                    h_d=sort(abs(h_d));
%                     p = 10^10*[P_T*abs(h_d(1)*h_d(2))^2, noise_power*(h_d'*h_d), -noise_power*abs(h_d(1))^2];
%                     r = roots(p);
%                     alpha=r(0<r&r<1); 
%                 h_d=sort(abs(h_d));
                final_rate=log2(1+P_T*h_d(1)^2/noise_power)/N1+final_rate; 
%                     final_rate=(1/N1)*(log2((abs(h_d(1))^2*P_T+noise_power)/(abs(h_d(1))^2*P_T*alpha+noise_power))+log2((abs(h_d(2))^2*alpha*P_T+noise_power)/noise_power))+final_rate;         
                elseif M_t>1
                    % I will write this later
                    [h_d, ~, ~]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                            x_i, y_i, z_i, alpha_d, alpha_r);      
                end 
                
            end
            final_rate=final_rate/N_iter;
        case{1} % with IRS
            
            if M_t==1
                for k=1:N_iter
                    %% channel gains
                    [h_d, g, h_r]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                              x_i, y_i, z_i, alpha_d, alpha_r);
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
                end
                    final_rate=final_rate/N_iter;
            elseif M_t>1
                % I will write this later
                    [h_d, G, h_r]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                              x_i, y_i, z_i, alpha_d, alpha_r);
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
end