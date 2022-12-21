% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% PotentialGame.m - This method returns the game table containing the utility function of the BS for each
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
    % interfernce        - Interfernce caused by IRS elements allocted to the other BS
%
% Outputs:
    % utility             - Utility function, 
    % rate                - Average acievable sumrate for each BS(|S1| * |S2| * |N|).
% ------------------------------------------------------------------------


function [utility, rate] = PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                            x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                            R, r, q_t, interference, epsilon)
   rate = zeros(ceil(1/q_t) + 1,ceil(1/q_t) + 1,2) ;
   utility = zeros(ceil(1/q_t) + 1,ceil(1/q_t) + 1,2) ;
   x0 = x ;
   y0 = y ;
      for j = 1 : N_iter
        t = 2 * pi * rand(N1,1) ;
        radius  =  R(1) * sqrt(rand(N1,1)) ;
        x(1 : N1) = x0(1 : N1)  +  radius.* cos(t) ;
        y(1 : N1) = y0(1 : N1)  +  radius.* sin(t)  ;
        t=2 * pi * rand(N2,1) ;
        radius = R(2) * sqrt(rand(N2,1)) ;
        x(N1 + 1 :  N1 + N2) = x0(N1 + 1 : N1 + N2)  +  radius.*cos(t) ;
        y(N1 + 1 :  N1 + N2) = y0(N1 + 1 : N1 + N2)  +  radius.*sin(t)  ;
        [h_d1, G1, h_r1] = ChannelGain(M_t, M_r, N1, N_i, x(1 : N1), y(1 : N1), z(1 : N1), x1, y1,...
                                                z1, x_i, y_i, z_i, alpha_d, alpha_r) ;
        [h_d2, G2, h_r2] = ChannelGain(M_t, M_r, N2, N_i, x(N1 + 1 : N1 + N2), y(N1 + 1 : N1 + N2),z(N1 + 1 : N1 + N2),...
                                    x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r) ;
        if M_t==1
          for s1 = 1 : ceil(1/q_t) + 1
               for s2 = 1 : ceil(1/q_t) + 2-s1
                N_i1 = floor((s1-1) * N_i * q_t) ;     % No. of elements allocated to BS1
                N_i2 = floor((s2-1) * N_i * q_t) ;     % No. of elements allocated to BS2
                [~, Z_cvx1] = Rate2(h_d1, G1(1 : N_i1), h_r1(1 : N_i1, : ), N1, N_i1, P_T, noise_power, N_iter) ;
                [~, Z_cvx2] = Rate2(h_d2, G2(1 : N_i2), h_r2(1 : N_i2, : ), N2, N_i2, P_T, noise_power, N_iter) ;
                Theta = diag([Z_cvx1 ; interference * Z_cvx2 ;zeros(N_i-N_i1-N_i2,1)]) ;
                channel_gain1 = sort(abs(h_d1 + h_r1' * Theta * G1)) ;
                [exact_rate1, ~] = SISONoma(epsilon, channel_gain1.^2 * P_T / noise_power) ;
%                 rate_lower_bound1 = N1 * log2(1 + P_T * channel_gain1(1)^2/noise_power)/N1 ;
%                 channel_gain1 = channel_gain.^2 * P_T / noise_power ; % normalized channel gain
                Theta = diag([Z_cvx2 ; interference * Z_cvx1 ;zeros(N_i-N_i1-N_i2,1)]) ;
                channel_gain2 = sort(abs(h_d2 + h_r2' * Theta * G2)) ;
                [exact_rate2, ~] = SISONoma(epsilon, channel_gain2.^2 * P_T / noise_power) ;
%                 rate_lower_bound2 = N2 * log2(1 + P_T * channel_gain2(1)^2/noise_power)/N2 ;
                
                rate(s1 , s2, 1) = rate(s1 , s2, 1) + ...
                                                    exact_rate1 ;              

                rate(s1 , s2, 2) = rate(s1 , s2, 2) + ...
                                                    exact_rate2 ;              
                                         
               end            
           end
        elseif M_t > 1
          for s1 = 1 : ceil(1/q_t) + 1
               for s2 = 1 : ceil(1/q_t) + 2-s1
                N_i1 = floor((s1-1) * N_i * q_t) ;     % No. of elements allocated to BS1
                N_i2 = floor((s2-1) * N_i * q_t) ;     % No. of elements allocated to BS2 
                [~ , exact_rate1, ~] = Rate3(h_d1, G1(1 : N_i1 , :), h_r1(1 : N_i1, : ), N1, N_i1, P_T, noise_power, epsilon) ;
                [~ , exact_rate2, ~] = Rate3(h_d2, G2(1 : N_i2 , :), h_r2(1 : N_i2, : ), N2, N_i2, P_T, noise_power, epsilon) ;          
                rate(s1 , s2, 1) = rate(s1 , s2, 1) + ...
                                                    exact_rate1 ;              
                rate(s1 , s2, 2) = rate(s1 , s2, 2) + ...
                                                    exact_rate2 ;              
               end
          end
        end
      end
        rate = rate/N_iter ;
end
