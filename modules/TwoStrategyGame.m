% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% TwoStrategyGame.m - This method returns the game table containing the
% utility function of the BS for each
%  strategy.
%                           Game Table:
%  ____________________________________________________________________
% |    |            with probability P(2)     |   1-P(2)               |
% |    |        W/O                           |     W                  |
% |___________________________________________|________________________|
% | W/O|    u(1,1,1),u(1,1,2)                 |    u(1,2,1),u(1,2,2)   |
% | W  |    u(2,1,1),u(2,1,2)                 |    u(2,2,1),u(2,2,2)   |
% |____|______________________________________|________________________|
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
    % interfernce        - Interfernce caused by IRS elements allocted to
    %                       the other BS
%
% Outputs:
    % utility             - Utility function, 
    % rate                - Average acievable rate for each BS's user(|S1| * |S2| * |N|).
% ------------------------------------------------------------------------


function [utility, rate] = TwoStrategyGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                            x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                            R, r, interference, epsilon)
    rate = zeros(2,2,2);
    utility = zeros(2,2,2) ;
    x0 = x ;
    y0 = y ;
    if interference==1
    if M_t==1
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
          for s1 = 1 : 2
               for s2 = 1 : 2
                N_i1 = floor((s1-1) * N_i) ;     % No. of elements allocated to BS1
                N_i2 = floor((s2-1) * N_i) ;     % No. of elements allocated to BS2
                if s1 == 2 & s2 == 2
                    N_i1 = N_i/2 ;
                    N_i2 = N_i/2 ;
                end
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
      end
    elseif M_t > 1
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
            for s1 = 1 : 2
                for s2 = 1 : 2
                    N_i1 = floor((s1-1) * N_i) ;     % No. of elements allocated to BS1
                    N_i2 = floor((s2-1) * N_i) ;     % No. of elements allocated to BS2
                    if s2 == 2 & s2 == 2
                        N_i1 = N_i/2 ;
                        N_i2 = N_i/2 ;
                    end
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
    elseif interference==0
    if M_t==1
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
          for s1 = 1 : 3
            N_i1 = floor((s1-1) * N_i)/2 ;     % No. of elements allocated to BS1
            [exact_rate1, Z_cvx1] = Rate2(h_d1, G1(1 : N_i1), h_r1(1 : N_i1, : ), N1, N_i1, P_T, noise_power, N_iter) ;
          if N1>2
            Theta = diag([Z_cvx1 ; zeros(N_i-N_i1,1)]) ;
            channel_gain1 = sort(abs(h_d1 + h_r1' * Theta * G1)) ;
            [exact_rate1, ~] = SISONoma(epsilon, channel_gain1.^2 * P_T / noise_power) ;
          end
            if s1 ==1
                rate(1 , :, 1) = rate(1 , :, 1) + exact_rate1 ;
            elseif s1==2
                rate(2 , 2, 1) = rate(2 , 2, 1) + exact_rate1 ;
            else 
                rate(2 , 1, 1) = rate(2 , 1, 1) + exact_rate1 ;
            end
          end
                for s2 = 1 : 3
                N_i2 = floor((s2-1) * N_i)/2 ;     % No. of elements allocated to BS2
                [exact_rate2, Z_cvx2] = Rate2(h_d2, G2(1 : N_i2), h_r2(1 : N_i2, : ), N2, N_i2, P_T, noise_power, N_iter) ;
                if N2>2
                Theta = diag([Z_cvx2 ; zeros(N_i-N_i2,1)]) ;
                channel_gain2 = sort(abs(h_d2 + h_r2' * Theta * G2)) ;
                [exact_rate2, ~] = SISONoma(epsilon, channel_gain2.^2 * P_T / noise_power) ;
                end
            if s2 ==1
                rate(: , 1, 2) = rate(: , 1, 2) + exact_rate2 ;
            elseif s2==2
                rate(2 , 2, 2) = rate(2 , 2, 2) + exact_rate2 ;
            else 
                rate(1 , 2, 2) = rate(1 , 2, 2) + exact_rate2 ;
            end
                end            
      end
    elseif M_t > 1
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
          for s1 = 1 : 3
            N_i1 = floor((s1-1) * N_i)/2 ;     % No. of elements allocated to BS1
            [~, Z_cvx1] = Rate2(h_d1, G1(1 : N_i1), h_r1(1 : N_i1, : ), N1, N_i1, P_T, noise_power, N_iter) ;
            Theta = diag([Z_cvx1 ; zeros(N_i-N_i1,1)]) ;
            channel_gain1 = sort(abs(h_d1 + h_r1' * Theta * G1)) ;
            [exact_rate1, ~] = SISONoma(epsilon, channel_gain1.^2 * P_T / noise_power) ;
            if s1 ==1
                rate(1 , :, 1) = rate(1 , :, 1) + exact_rate1 ;
            elseif s1==2
                rate(2 , 2, 1) = rate(2 , 2, 1) + exact_rate1 ;
            else 
                rate(2 , 1, 1) = rate(2 , 1, 1) + exact_rate1 ;
            end
          end
        for s2 = 1 : 3
                N_i2 = floor((s2-1) * N_i)/2 ;     % No. of elements allocated to BS2
                [~, Z_cvx2] = Rate2(h_d2, G2(1 : N_i2), h_r2(1 : N_i2, : ), N2, N_i2, P_T, noise_power, N_iter) ;
%                 rate_lower_bound1 = N1 * log2(1 + P_T * channel_gain1(1)^2/noise_power)/N1 ;
%                 channel_gain1 = channel_gain.^2 * P_T / noise_power ; % normalized channel gain
                Theta = diag([Z_cvx2 ; zeros(N_i-N_i2,1)]) ;
                channel_gain2 = sort(abs(h_d2 + h_r2' * Theta * G2)) ;
                [exact_rate2, ~] = SISONoma(epsilon, channel_gain2.^2 * P_T / noise_power) ;
%                 rate_lower_bound2 = N2 * log2(1 + P_T * channel_gain2(1)^2/noise_power)/N2 ;
            if s2 ==1
                rate(: , 1, 2) = rate(: , 1, 2) + exact_rate2 ;
            elseif s2==2
                rate(2 , 2, 2) = rate(2 , 2, 2) + exact_rate2 ;
            else 
                rate(1 , 2, 2) = rate(1 , 2, 2) + exact_rate2 ;
            end
        end            
      end
end
    end
rate = rate/N_iter ;
