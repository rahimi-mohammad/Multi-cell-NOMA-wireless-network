% ------------------------------------------------------------------------
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method calculates optimal acheivable rate and IRS phase shift
% of the BS's users
% for a given beamforming
% Inputs:
    % W                  - Beamforming  [w1 , w2 , ... , w_N1 ]
    % epsilon            - 
    % h_d                - BS to users channel gain(M_t * M_r * N1) , 
    % G                  - BS to IRS channel gain , 
    % h_r                - IRS to users channel gain(N * N1).    
% Outputs:
    % rate               - Acheivable Rate of each user(the exact rate).
    % b_star             - Users optimal power allocation.
% ------------------------------------------------------------------------
function [rate, Z_cvx , feasible] = MISONoma1(W , epsilon, h_d , G , h_r , noise_power)
feasible = false ;    
q_min = 0 ;
    q_max = 10000 ;
    N1 = size(h_d , 2) ; % No. users
    N = size(h_r , 1) ; % No. IRS elements
    M_t = size(h_d , 1) ; % No. transmitter antennas
    %% problem coefficients
    R = zeros( N + 1 , N + 1 , N1 , N1);
    S = zeros( N + 1 , N + 1, N1 );
    J = zeros( N , M_t , N1) ;  %  N * M_t
    for k = 1 : N1
        J( : , : , k ) = diag(h_r(: , k)') * G ;     % N * M_t
        S( : , : ,k) =1e10 * [J(: , : , k) * J(: , : , k)' , J(: , : , k) * h_d(  : , k) ;...
            h_d( : , k )' * J(: , : , k)'   ,0] ; 
        for t = 1 : N1 
            l = diag(h_r(: , k)') * G * W(: , t) ;
            v = h_d( : , k)' * W(: , t) ;
            R( : , : ,  k ,  t ) =1e10 * [ l * l'  ,  l * v'  ; l' * v , 0 ];
        end
    end
%% problem solver
    X_star = zeros(N + 1,N + 1) ;
    while q_max - q_min >= epsilon
        q = (q_min + q_max)/2 ;
        cvx_begin sdp quiet
%             dual variable y
            variable X(N + 1,N + 1) complex hermitian  ;
            minimize -q 
            subject to 
            for k = 1 : N1
                for t = 1 : k
                    A = sum(R( : , : , k , t + 1 : N1) , 4 ) ;
                    
                   1e10 * q * noise_power + q * trace(  A * X ) + ...
                   1e10 * q * sum(abs(h_d( : , k)' * W(: , t + 1 : N1 )).^2)  - ...
                   trace(R( : , : , k , t) * X ) -...
                   1e10 * abs(h_d(:, k)' * W(: , t))^2 <= 0
                end
            end
            for k = 1 : N1 - 1
                1 *real(trace(S( : , : ,k) * X)) +...
                    1e10 *real(h_d(: , k)' * h_d(: , k)) - ...
                    1 *real(trace(S( : , : ,k + 1) * X)) -...
                    1e10 *real(h_d(: , k + 1)' * h_d(: , k + 1))  <= 0
            end
            X == hermitian_semidefinite(N + 1) ; 
            diag(X) == 1 ; 
        cvx_end
        if cvx_status( 1 : 6 ) == 'Solved'
%             disp('It is feasible')
            feasible = true ;
            X_star = X ;
            q_min = q ;
        else
            q_max = q ;
        end
    end
    X = X_star ;
    Z_cvx = GR(X , 1) ;
    rate = log2(1 + q) ;