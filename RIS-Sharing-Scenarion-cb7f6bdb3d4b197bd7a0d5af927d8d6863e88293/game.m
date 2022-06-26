% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% game.m - This method returns the game table containing the utility function of the BS for each
%  strategy.
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
%
%
%
%
%
%
%
% Outputs:
    % utility                  - Utility function, 
    % rate               - Average acievable rate for each BS.
% ------------------------------------------------------------------------


function [utility, rate]=game(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter, R, r)
   rate=zeros(2,2,2);
   x0= x;
   y0= y;
   fig1=figure;
   for j=1:N_iter
        s=1;                            % s=1 : with IRS
        t=2*pi*rand(N1,1);
        radius = R(1)*sqrt(rand(N1,1));
        x(1:N1) = x0(1:N1) + radius.*cos(t);
        y(1:N1) = y0(1:N1) + radius.*sin(t) ;
        t=2*pi*rand(N2,1);
        radius = R(2)*sqrt(rand(N2,1));
        x(N1+1: N1+N2) = x0(N1+1:N1+N2) + radius.*cos(t);
        y(N1+1: N1+N2) = y0(N1+1:N1+N2) + radius.*sin(t) ;
        plot(x(1),y(1),'ro','linewidth',1)
        hold on
        plot(x(2),y(2),'bo','linewidth',1)
        plot(x(3),y(3),'ro','linewidth',1)
        plot(x(4),y(4),'bo','linewidth',1)
        plot(x(5),y(5),'ro','linewidth',1)
        rate(2, 2, 1) = rate(2, 2, 1) + Rate(s, M_t, M_r, N1, N_i/2, P_T, x(1: N1), y(1: N1),...
                                            z(1: N1), x1, y1, z1, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
        rate(2, 2, 2) = rate(2, 2, 2) + Rate(s,  M_t, M_r, N2, N_i/2, P_T, x(N1+1: N1+N2), y(N1+1: N1+N2),...
                                            z(N1+1: N1+N2), x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
        rate(1, 2, 2) = rate(1, 2, 2) + Rate(s, M_t, M_r, N2, N_i, P_T, x(N1+1: N1+N2), y(N1+1: N1+N2),...
            z(N1+1: N1+N2), x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);    
        
       
        rate(2, 1, 1) = rate(2, 1, 1) + Rate(s, M_t, M_r, N1, N_i, P_T, x(1: N1), y(1: N1),....
                                            z(1: N1), x1, y1, z1, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
        s=0;                           % s=0 : without IRS
        rate(1 , 1, 1) = rate(1, 1, 1) + Rate(s, M_t, M_r, N1, N_i/2, P_T, x(1:N1), y(1:N1), ...
                                            z(1:N1), x1, y1, z1, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
        rate(1, 1, 2) = rate(1, 1, 2) + Rate(s, M_t, M_r, N2, N_i/2, P_T, x(N1+1: N1+N2), y(N1+1: N1+N2), ...
                                            z(N1+1: N1+N2), x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, 1);
        rate(1, 2, 1) = rate(1, 1, 1);
        rate(2, 1, 2) = rate(1, 1, 2);
   end
   close(fig1)
        rate=rate/N_iter;
        utility(1, 1, 1)=rate(1, 1, 1);
        utility(1, 1, 2)=rate(1, 1, 2);
        utility(1, 2, 1)=rate(1, 2, 1);    
        utility(2, 1, 2)=rate(2, 1, 2);

        utility(1, 2, 2)=rate(1, 2, 2)-r*N_i;
        utility(2, 1, 1)=rate(2, 1, 1)-r*N_i;

        utility(2, 2, 1)=utility(2, 1, 1)-r*N_i/2;
        utility(2, 2, 2)=utility(1, 2, 2)-r*N_i/2;        
end
