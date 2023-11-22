% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% StackelbergGame.m - This method returns the game table containing the utility function of the BS for each
%  strategy.
%                                     _____________
%                                    |Leader : IRS |
%                                    |_____________|
%                                   /               \
%                                  /                 \
%                     _____________                      _____________
%       Followers:   |     BS1     |                    |      BS2    |
%                    |_____________|                    |_____________|
%                                             
% Inputs:
    % N1                 - No.  users of BS1,
    % N_i                - No.  IRS elements,
    % P_T                - BS power,
    % x,y,z              - Users location,
    % x1,y1,z1           - BS location,   
    % x_i,y_i,z_i        - IRS location,      
    % alpha_d            - BS to user pathloss,
    % alpha_r            - BS to IRS pathloss,
    % noise_power        - Noise power,
    % N_iter             - No. iteraions,
    % R                  - User area center radius,
    % q_t                - Quantization step-size,
    % interfernce        - Interfernce caused by IRS elements allocted to
    %                      the other BS,
    % StepSize           - Step size in algorithm.
% Outputs:
    % r                   - Price of IRS per element, 
    % rate                - Average acievable sumrate for each BS.
% ------------------------------------------------------------------------


function [r, rate, N_i_star] = StackelbergGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                            x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                            R, q_t, interference, StepSize, r_0)
   r = r_0 ;
   N_i_0 = zeros(2, 1) ;
   N_i_star = zeros(2, 1) ;
   epsilon  = 0.002 ;
%    utility = zeros(ceil(1/q_t) + 1,ceil(1/q_t) + 1,2) ;
   EPG_utility = zeros(ceil(1/q_t)+1,ceil(1/q_t)+1,2) ; 
    s1_star = 0 ;
    s2_star = 0 ;
   [ ~ , rate] = PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                            x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                            R, r, q_t, interference) ;
   N_i1 = 0 ;
   N_i2 = 0 ;
   while 'True'
       for s1 = 1:ceil(1/q_t)+1
           for s2 = 1:ceil(1/q_t)+1
               N_i1 = floor((s1-1) * N_i * q_t) ;     % No. of elements allocated to BS1
               N_i2 = floor((s2-1) * N_i * q_t) ;     % No. of elements allocated to BS2
               EPG_utility(s1, s2, 1) = rate(s1, s2, 1) - N_i1 * r ; 
               EPG_utility(s1, s2, 2) = rate(s1, s2, 2) - N_i2 * r ; 
           end            
       end
        F = sum(EPG_utility,3) ; 
        maximum = max(max(F)) ; 
        [s1_star,s2_star] = find(F == maximum) ;
       N_i_0 = N_i_star ;
       N_i_star = [floor((s1_star-1) * N_i * q_t) ; floor((s2_star-1) * N_i * q_t) ]      % No. of elements allocated to BS2
       N_i1 = [ N_i1 ; N_i_star(1) ] ;
       N_i2 = [ N_i2 ; N_i_star(2) ] ;
%        r_0 = r ;
%        r = r - StepSize * sign( N_i/2 - sum(N_i_star ,"all"))
%        r = r + StepSize * sign(1 + sum(N_i_0 ,"all") - sum(N_i_star ,"all"))
       r = r + StepSize * 1 ;
%        if sign(1 + sum(N_i_0 ,"all") - sum(N_i_star ,"all")) < 0
%             StepSize= 0.5 * StepSize ;
%        end
%        abs(r * sum(N_i_star ,"all") - r_0 * sum(N_i_0 ,"all"))
%        if abs(r * sum(N_i_star ,"all") - r_0 * sum(N_i_0 ,"all")) < epsilon
%            break ;
%        end
        if sum( N_i_star , "all" ) == 0
            break ;
        end
    end
plot( r_0 : StepSize : r , N_i1(2 : end) )
hold on
plot( r_0 : StepSize : r , N_i1(2 : end) )
legend('N_i1' , 'N_i2')
%    rate = rate(s1_star,s2_star,:) ; 
end
