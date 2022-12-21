% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method calculates optimal acheivable rate and BSs' beamforming
% for a given IRS phase shift.
% Inputs:
    % epsilon            - epsilon for Bisection method
    % h_d                - BS to users channel gain(M_t * M_r * N1) , 
    % noise_power        - Noise power,
    % P_T                - BS power.
% Outputs:
    % rate               - Acheivable Rate of each user(the exact rate).
    % b_star             - Users optimal power allocation.
% ------------------------------------------------------------------------

function [rate, W] = MISONoma2( epsilon, h_d , noise_power , P_T)
    q_min = 0 ;
    q_max = 10000 ;
    M_t = size(h_d , 1) ; % No. transmitter antennas
    N1 = size(h_d , 2) ;  % No. users
    W = zeros( M_t , N1) ;
    h = h_d ;
    %% problem coefficients
    H = zeros(  M_t , M_t , N1);
    for k = 1 : N1
        H( : , : , k) = 1e10 * h(: , k) * h(: , k)' ;
    end
%% problem solver
    while q_max - q_min >= epsilon

        q = (q_min + q_max)/2 ;
        cvx_begin sdp quiet
            variable X( M_t  , M_t , N1 ) complex hermitian  ;
            minimize -q 
            subject to 
            
            for k = 1 : N1 - 1
                X( : , : , k) == hermitian_semidefinite(M_t) ; 
                diag(X( : , : , k)) >= 0
                for t = 1 : k 
                    1e10 * q * noise_power +   1 * q * real(trace( H( : , : , k )*...
                        sum( X(: , : , t + 1 : N1) , 3 )  ) - ...
                            1 * trace(H( : , : , k ) * X(: , : , t))) <= 0
                end
            end
            k = N1 ;
            X( : , : , k) == hermitian_semidefinite(M_t) ; 
            diag(X( : , : , k)) >= 0
            for t = 1 : k - 1 
                1e10 * q * noise_power +...
                1 * q * real(trace(  H( : , : , k )* sum( X(: , : , t + 1 : N1) , 3 )  )) - ...
                1 * real(trace(H( : , : , k ) * X(: , : , t))) <= 0
            end
            1e10 * q * noise_power - 1 * real(trace(H( : , : , k ) * X(: , : , k))) <= 0
            trace(sum(X, 3)) <= P_T
            
        cvx_end
        
        if cvx_status( 1 : 6 ) == 'Solved'
%             disp('It is feasible')
            X_star = X ;
            q_min = q ;
        else
            q_max = q ;
        end
    end
    for k = 1 : N1
            W( : , k) = GR(X_star( : , :,k) , 2) ;
    end
    q = q_min ;
    rate = log2(1 + q) ;