function u=Utility(s,N1,N_i,P_T,x,y,z,x1,y1,z1,x_i,y_i,z_i,alpha_d,alpha_r,noise_power)
    %% parameters
    M_t=1;                          % No. of transmitter antennas
    M_r=1;                          % No. of receiver antennas
    N_iter=100;                     % No. iteraions
    channel_gain=zeros(N1,1);
    final_rate=0;
    N=N_i;
        m=repmat('.', 1, 11);
    switch s
        case{0} % without RIS
            for k=1:N_iter
                m(1+floor(10*k/N_iter))='|';
                disp(m)
                %% channel gains
                h_d=(10^0.2)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x-x1).^2+(y-y1).^2+(z-z1).^2 )).^alpha_d);
                %% without IRS
                h_d=sort(abs(h_d));
                p = 10^10*[P_T*abs(h_d(1)*h_d(2))^2, noise_power*(h_d'*h_d), -noise_power*abs(h_d(1))^2];
                r = roots(p);
                alpha=r(0<r&r<1); 
                final_rate=(1/N1)*(log2((abs(h_d(1))^2*P_T+noise_power)/(abs(h_d(1))^2*P_T*alpha+noise_power))+log2((abs(h_d(2))^2*alpha*P_T+noise_power)/noise_power))+final_rate;         
                
            end
            final_rate=final_rate/N_iter;
        case{1}
            for k=1:N_iter
                m(1+floor(10*k/N_iter))='|';
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
                    maximize min(10^13*abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),10^13*abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
                    subject to 
                    diag(X) == 1;
                    X==hermitian_semidefinite(N+1);
                cvx_end
                Z_cvx=GR(X);
%                 Q=10^(-13)*Q;
                %% channel gains
                channel_gain(1)=abs(h_d(1)+h_r(:,1)'*diag(Z_cvx)*g);
                channel_gain(2)=abs(h_d(2)+h_r(:,2)'*diag(Z_cvx)*g);
                channel_gain=sort(abs(channel_gain));
                p = 10^12*[P_T*abs(channel_gain(1)*channel_gain(2))^2, noise_power*(channel_gain'*channel_gain), -noise_power*abs(channel_gain(1))^2];
                r = roots(p);
                alpha=r(0<r&r<1); 
                rate_lower_bound=(1/N1)*(log2((abs(channel_gain(1))^2*P_T+noise_power)/(abs(channel_gain(1))^2*P_T*alpha+noise_power))+log2((abs(channel_gain(2))^2*alpha*P_T+noise_power)/noise_power));
                final_rate=final_rate+rate_lower_bound;
            end
                final_rate=final_rate/N_iter;
    end
    u=final_rate;
end