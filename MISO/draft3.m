%% Test Beamforming
% 
%   
%% 
% clc

h_d = h_d1 ;
G = G1 ;
h_r = h_r1 ;
W = 1 * randn(M_t , N1)/sqrt(4 * N1 * M_t/ P_T) ; 
h = zeros(M_t , N1) ; % Cascade channel
        q_min = 1 ;
    q_max = 10000 ;
    N1 = size(h_d , 2) ; % No. users
    N = size(h_r , 1) ; % No. IRS elements
    M_t = size(h_d , 1) ; % No. transmitter antennas
    %% problem coefficients
    R = zeros( N + 1 , N + 1 , N1 , N1);
    S = zeros( N + 1 , N + 1, N1 );
    J = zeros( N , M_t , N1) ;  %  N * M_t
    for k = 1 : N1
        J( : , : , k ) = diag(h_r(: , k)') * G ; % N * M_t
        S( : , : ,k) =1e10 * [J(: , : , k) * J(: , : , k)' , J(: , : , k) * h_d(  : , k) ;...
            h_d( : , k )' * J(: , : , k)'   ,0] ; 
        for t = 1 : k 
            l = diag(h_r(: , k)') * G * W(: , t) ;
            v = h_d( : , k)' * W(: , t) ;
            R( : , : ,  k ,  t ) =1e10 * [ l * l'  ,  l * v'  ; l' * v , 0 ];
        end
    end
%% problem solver
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
                    trace(R( : , : , k , t) * X) -...
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
            X_star = X ;
            q_min = q ;
        else
            q_max = q ;
        end
    end
    X = X_star ;
    Z_cvx = GR(X , 1) ;
    rate = log2(q) ;
    
%% %%% test
for k = 1 : N1
    h( : , k ) = (h_r(: , k)' * diag(Z_cvx) * G + h_d(: , k)')' ;
end
r1 = log2(1 + (abs(h(:, 1)' * W(: , 1))^2)/(noise_power + abs(h(:, 1)' * W(: , 2))^2 )) ;
r21 = log2(1 + (abs(h(:, 2)' * W(: , 1))^2)/(noise_power + abs(h(:, 2)' * W(: , 2))^2 )) ; 
r22 = log2(1 + (abs(h(:, 2)' * W(: , 2))^2)/noise_power) ; 
disp(r1)
disp(r21)
disp(r22)  
    
    
%% %%%%


   

%                 rate0 = 0 ;
%             rate2 = 0.2 ;
%     while abs(rate2 - rate0) > 0.1
%             rate0 = rate2 ;
%             [rate1, Z_cvx] = MISONoma1(W , epsilon, h_d , G , h_r , noise_power) ;  
%             rate1
%                 for k = 1 : N1
%                     h( : , k ) = (h_r(: , k)' * diag(Z_cvx) * G + h_d(: , k)')' ;
%                 end
%             %% Solve P1.4
%             [rate2, W] = MISONoma2(epsilon, h , noise_power , P_T) ;
%             disp('after beamforming')
%             rate2
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
%             
% r1 = log2(1 + (abs(h(:, 1)' * W(: , 1))^2)/(noise_power + abs(h(:, 1)' * W(: , 2))^2 )) ;
% r21 = log2(1 + (abs(h(:, 2)' * W(: , 1))^2)/(noise_power + abs(h(:, 2)' * W(: , 2))^2 )) ; 
% r22 = log2(1 + (abs(h(:, 2)' * W(: , 2))^2)/noise_power) ; 
% disp(r1)
% disp(r21)
% disp(r22)  
%     
%    
%     
%     
%     
%     
%     
%     
%     
%     
%     
    
    
