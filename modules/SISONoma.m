% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method calculates optimal aceivable rate of the BS's users
% for a given channel gain
% Inputs:
    % epsilon            - 
    % channel_gain       - Normalized BS to users channel gain,
% Outputs:
    % rate               - Acheivable Rate of each user(the exact rate).
    % b_star             - Users optimal power allocation.
% ------------------------------------------------------------------------
function [rate, b_star] = SISONoma(epsilon, channel_gain)
%% calculating optimal rate
q_min = 1 ;
q_max = 1000 ;
N = size(channel_gain , 1) ;
b_star = zeros( N+1, 1) ;
    if N==2
        q_star = (channel_gain(1) - channel_gain(2) +...
            sqrt( (channel_gain(1) - channel_gain(2))^2 +...
            4*channel_gain(1)*(channel_gain(2) +...
            channel_gain(1)*channel_gain(2)) ))/2/channel_gain(1) ;
        rate = log2(q_star)  ;
        b_star(2 , 1) = ( q_star - 1 )/channel_gain(2) ;
        b_star(1 , 1) = 1- b_star(2 , 1) ;
    else
        while q_max - q_min >= epsilon
            q = (q_min + q_max)/2 ;
            cvx_begin quiet
                dual variable y
                variable b(N + 1)   ;
                minimize -q 
                subject to 
                b(1) == 1 ;
                b(N + 1) == 0 ;
                for m = 1 : N 
                    b(m + 1) - b(m) <= 0
                    -b(m) <= 0
                    -(1 + b(m) * channel_gain(m)) + q + b(m + 1) * channel_gain(m) * q <= 0
                end
            cvx_end
            if cvx_status( 1 : 6 ) == 'Infeas'
                q_max = q ;
            else
                q_min = q ;
                b_star = b ;
            end
        end
        rate = log2(q) ;
    end
end



