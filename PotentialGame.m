% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% game.m - This method returns the game table containing the utility function of the BS for each
%  strategy.
%                           Game Table:
%  ______________________________________________________________________
% |    |  0                   |        ...      |     N(tiles)           |
% |____|______________________|_______________ _|________________________|
% |  0 |    u(1,1,1),u(1,1,2) |                 | u(1,N+1,1),u(1,N+1,2)  |
% |____|______________________|_________________|________________________|
% | .  |       .              | .               |       .                |
% | .  |       .              |    .            |       .                |
% | .  |       .              |       .         |       .                |
% |____|______________________|_________________|________________________|
% | N  |u(N+1,1,1),u(N+1,1,2) |                 | u(N+1,N+1,1),u(2,2,2)  |
% |____|_____________________ |_________________|________________________|
%-------------------------------------------------------------------------
% Inputs:
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
    % R                  - User area center radius
    % q_t                - Quantization step-size

%
% Outputs:
    % utility             - Utility function, 
    % rate                - Average acievable sumrate for each BS.
% ------------------------------------------------------------------------


function [utility, rate]=PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter, R, r, q_t)
   rate=zeros(ceil(1/q_t)+1,ceil(1/q_t)+1,2);
   utility=zeros(ceil(1/q_t)+1,ceil(1/q_t)+1,2);
   x0= x;
   y0= y;
   
   for j=1:N_iter
%         s=1;                            % s=1 : with IRS
        t=2*pi*rand(N1,1);
        radius = R(1)*sqrt(rand(N1,1));
        x(1:N1) = x0(1:N1) + radius.*cos(t);
        y(1:N1) = y0(1:N1) + radius.*sin(t) ;
        t=2*pi*rand(N2,1);
        radius = R(2)*sqrt(rand(N2,1));
        x(N1+1: N1+N2) = x0(N1+1:N1+N2) + radius.*cos(t);
        y(N1+1: N1+N2) = y0(N1+1:N1+N2) + radius.*sin(t) ;
%         plot(x(1),y(1),'ro','linewidth',1)
%         hold on
%         plot(x(2),y(2),'bo','linewidth',1)
%         plot(x(3),y(3),'ro','linewidth',1)
%         plot(x(4),y(4),'bo','linewidth',1)
%         plot(x(5),y(5),'ro','linewidth',1)
% 
%     temp=Rate(M_t, M_r, N1, transpose(0:N_i*q_t:N_i), P_T, x(1:N1), y(1:N1), ...
%         z(1:N1), x1, y1, z1, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
% 
%         rate(s1 , 1:ceil(1/q_t)+2-s1, 1) + Rate(M_t, M_r, N1, transpose(0:N_i*q_t:N_i), P_T, x(1:N1), y(1:N1), ...
%                                                 z(1:N1), x1, y1, z1, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);

        for s1=1:ceil(1/q_t)+1
            N_i1=floor((s1-1)*N_i*q_t);     % No. of elements allocated to BS1
            rate(s1 , 1:ceil(1/q_t)+2-s1, 1) = rate(s1 , 1:ceil(1/q_t)+2-s1, 1) +...
                                            Rate( M_t, M_r, N1, N_i1,...
                                            P_T, x(1:N1), y(1:N1), z(1:N1), x1, y1,...
                                            z1, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
           for s2=1:ceil(1/q_t)+1-s1
                N_i2=floor((s2-1)*N_i*q_t);    % No. of elements allocated to BS2
                rate(s1 , s2, 2) = rate(s1 , s2, 2) + ...
                                    Rate(M_t, M_r, N2, N_i2, P_T,...
                                    x(N1+1: N1+N2), y(N1+1: N1+N2),z(N1+1: N1+N2),...
                                    x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);                          
           end            
        end
        
   end
   
        rate=rate/N_iter;
        utility=utility/N_iter;
        rate(:, :, 1)=N1*rate(:, :, 1);
        rate(:, :, 2)=N2*rate(:, :, 2);

        for s1=1:ceil(1/q_t)+1
           for s2=1:ceil(1/q_t)+1
                N_i1=floor((s1-1)*N_i*q_t);    
                N_i2=floor((s2-1)*N_i*q_t);    % No. of elements allocated to BS2
               utility(s1 , s2, 1)=rate(s1 , s2, 1)-r*N_i1;
               utility(s1 , s2, 2)=rate(s1 , s2, 2)-r*N_i2;
           end            
        end
        
%     for s1=1:N_i*q_t+1
%            for s2=1:N_i*q_t-s1+1
%                N_i1=floor((s1-1)/q_t);     % No. of elements allocated to BS1
%                N_i2=floor((s2-1)/q_t);    % No. of elements allocated to BS2
%                utility(s1 , s2, 1)=rate(s1 , s2, 1)-r*N_i1;
%                utility(s1 , s2, 2)=rate(s1 , s2, 2)-r*N_i2;           
%            end            
%     end
        
end
