function [rate_lower_bound1,rate_lower_bound2]=phase_shift(s,N1,N,x1,y1,z1,x2,y2,z2,x_i,y_i,z_i,alpha_d,alpha_r,M_t,M_r,x,y,z,P_T,noise_power) 
            h_d1=(10^0.2)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(1:2)-x1).^2+(y(1:2)-y1).^2+(z(1:2)-z1).^2 )).^alpha_d);
            h_d2=(10^0.2)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(3:4)-x2).^2+(y(3:4)-y2).^2+(z(3:4)-z2).^2 )).^alpha_d);
         if s==1
            N=floor(N/2);
              
              %% channel gains1
            g1=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^alpha_r);
            h_r1=zeros(N,N1);
            for i=1:N1
                h_r1(:,i)=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha_r);
            end
              %% channel gains2
            g2=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x2).^2+(y_i-y2).^2+(z_i-z2).^2 )).^alpha_r);
            h_r2=zeros(N,N1);
            for i=1:N1
                h_r2(:,i)=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i+2)-x_i).^2+(y(i+2)-y_i).^2+(z(i+2)-z_i).^2 )).^alpha_r);
            end
            %% problem coefficients
            phi1=zeros(N,N1);
            Q1=zeros(N+1,N+1,N1);
            for i=1:N1
                phi1(:,i)=diag(h_r1(:,i)')*g1;
                Q1(:,:,i)=[phi1(:,i)*transpose(conj(phi1(:,i))) conj(h_d1(i))*phi1(:,i);h_d1(i)*transpose(conj(phi1(:,i))) 0];
            end
            Q1=10^10*Q1;
            %% problem solver1
            cvx_begin sdp
                variable X(N+1,N+1) complex hermitian  ;
                maximize min(10^10*abs(h_d1(1))^2+real(trace(Q1(:,:,1)*X)),abs(10^10*h_d1(2))^2+real(trace(Q1(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
                subject to 
                diag(X) == 1;
                X==hermitian_semidefinite(N+1);
            cvx_end
            Z_cvx1=GR(X);
            %% problem coefficients
            phi2=zeros(N,N1);
            Q2=zeros(N+1,N+1,N1);
            for i=1:N1
                phi2(:,i)=diag(h_r2(:,i)')*g2;
                Q2(:,:,i)=[phi2(:,i)*transpose(conj(phi2(:,i))) conj(h_d2(i))*phi2(:,i);h_d2(i)*transpose(conj(phi2(:,i))) 0];
            end
            Q2=10^10*Q2;
            %% problem solver2

            cvx_begin sdp
                variable X(N+1,N+1) complex hermitian  ;
                maximize min(10^10*abs(h_d2(1))^2+real(trace(Q2(:,:,1)*X)),abs(10^10*h_d2(2))^2+real(trace(Q2(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
                subject to 
                diag(X) == 1;
                X==hermitian_semidefinite(N+1);
            cvx_end
            Z_cvx2=GR(X); 
            

            channel_gain1(1)=abs(h_d1(1)+h_r1(:,1)'*diag(Z_cvx1)*g1);
            channel_gain1(2)=abs(h_d1(2)+h_r1(:,2)'*diag(Z_cvx1)*g1);
            
            channel_gain2(1)=abs(h_d2(1)+h_r2(:,1)'*diag(Z_cvx2)*g2);
            channel_gain2(2)=abs(h_d2(2)+h_r2(:,2)'*diag(Z_cvx2)*g2);   
            
            rate_lower_bound1=(1/N1)*(log2(1+min(channel_gain1.^2)*P_T/noise_power));
            rate_lower_bound2=(1/N1)*(log2(1+min(channel_gain2.^2)*P_T/noise_power));
         end
        
         
         if s==2
             %% channel gains1
            g1=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^alpha_r);
            h_r1=zeros(N,N1);
            for i=1:N1
                h_r1(:,i)=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha_r);
            end
            %% problem coefficients
            phi1=zeros(N,N1);
            Q1=zeros(N+1,N+1,N1);
            for i=1:N1
                phi1(:,i)=diag(h_r1(:,i)')*g1;
                Q1(:,:,i)=[phi1(:,i)*transpose(conj(phi1(:,i))) conj(h_d1(i))*phi1(:,i);h_d1(i)*transpose(conj(phi1(:,i))) 0];
            end
            Q1=10^10*Q1;
            %% problem solver1
            cvx_begin sdp
                variable X(N+1,N+1) complex hermitian  ;
                maximize min(10^10*abs(h_d1(1))^2+real(trace(Q1(:,:,1)*X)),abs(10^10*h_d1(2))^2+real(trace(Q1(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
                subject to 
                diag(X) == 1;
                X==hermitian_semidefinite(N+1);
            cvx_end
            Z_cvx1=GR(X);             
             
            channel_gain1(1)=abs(h_d1(1)+h_r1(:,1)'*diag(Z_cvx1)*g1);
            channel_gain1(2)=abs(h_d1(2)+h_r1(:,2)'*diag(Z_cvx1)*g1);
            
            channel_gain2(1)=abs(h_d2(1)+h_r2(:,1)'*diag(Z_cvx2)*g2);
            channel_gain2(2)=abs(h_d2(2)+h_r2(:,2)'*diag(Z_cvx2)*g2);   
            
            rate_lower_bound1=(1/N1)*(log2(1+min(channel_gain1.^2)*P_T/noise_power));
            rate_lower_bound2=(1/N1)*(log2(1+min(channel_gain2.^2)*P_T/noise_power));
             
             
         end
    
end