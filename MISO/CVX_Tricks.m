%% CVX Tricks
% Remember that the matrix multiplied by you variable can not be too small
% in quantity of its elements. look at these two examples. Here matrix H
% has smal elements
%% 
N1 = 2 ;                                 % No. first Basestations users
R1 = 10 ;                                % BS1 user area center radius
R2 = 10 ;                                % BS2 user area center radius    
N = 5 ;                                % No.  IRS elements
P_T = (10^(40/10)) * 1e-3 ;                 % BS power
[x1 , y1 , z1] = deal(0 , 0 , 0) ;                % BS1 location
[x2 , y2 , z2] = deal(100 , 10 , 0) ;             % BS2 location
[x_i , y_i , z_i] = deal(50 , 30 , 10) ;          % IRS location
M_t = 2 ;                                 % No. of transmitter antennas
M_r = 1 ;
d_IB = 100 ; 
alpha_d = 3.76 ; 
alpha_r = 2 ; 
noise_power = (10^(-114/10)) ;            % -169dbm/Hz
epsilon = 0.05 ;
%% users center location
x0 = zeros(N1+N2 , 1) ; 
y0 = zeros(N1+N2 , 1) ; 
z = 1.5 * ones(N1+N2 , 1) ; 
x0(1 : N1 , 1) = 35 ; 
x0(N1+1 : N2+N1) = 55 ; 
y0(1 : N1 , 1) = 25 ; 
y0(N1+1 : N2+N1) = 10 ;
R = [10 ; 10] ;
[h_d, G, h_r] = ChannelGain(M_t, M_r, N1, N, x(1 : N1 , 1), y(1 : N1 , 1), z(1 : N1), x1, y1,...
                                    z1, x_i, y_i, z_i, alpha_d, alpha_r) ;
Z_cvx = exp(1j * rand( N , 1) * 2 * pi) ;
for k = 1 : N1
    h( : , k ) = (h_r(: , k)' * diag(Z_cvx) * G + h_d(: , k)')' ;
end

M_t = size(h_d , 1) ; % No. transmitter antennas
N1 = size(h_d , 2) ;  % No. users
W = zeros( M_t , N1) ;
H = zeros(  M_t , M_t , N1);
    
%% Example 1
q_min = 1 ;
q_max = 10000 ;
    for k = 1 : N1
        H( : , : , k) = 1e10 * h(: , k) * h(: , k)' ;
    end
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
    rate = log2(q) ;
X = X_star ;
H = H * 1e-10 ;
% Evaluation of example 1 without GR
r1 = log2(1 + (1 * trace(H( : , : , 1 ) * X(: , : , 1))) /...
    (noise_power+ real(trace( H( : , : , 1 )*...
    sum( X(: , : , 2) , 3 )  ) ))) ;
r21 = log2(1 + (1 * trace(H( : , : , 2 ) * X(: , : , 1))) /...
    (noise_power+ real(trace( H( : , : , 2 )*...
    sum( X(: , : , 2) , 3 )  ) ))) ;
r22 = log2(1 + (1 * trace(H( : , : , 2 ) * X(: , : , 2))) / noise_power) ;
disp(r1)
disp(r21)
disp(r22)

% Evaluation of example 1 after GR
r1 = log2(1 + (abs(h(:, 1)' * W(: , 1))^2)/(noise_power + abs(h(:, 1)' * W(: , 2))^2 )) ;
r21 = log2(1 + (abs(h(:, 2)' * W(: , 1))^2)/(noise_power + abs(h(:, 2)' * W(: , 2))^2 )) ; 
r22 = log2(1 + (abs(h(:, 2)' * W(: , 2))^2)/noise_power) ; 
disp(r1)
disp(r21)
disp(r22)




%  %% Example 2
%      q_min = 1 ;
%     q_max = 10000 ;
%     for k = 1 : N1
%         H( : , : , k) = 1 * h(: , k) * h(: , k)' ;
%     end
%     while q_max - q_min >= epsilon
% 
%         q = (q_min + q_max)/2 ;
%         cvx_begin sdp quiet
%             variable X( M_t  , M_t , N1 ) complex hermitian  ;
%             minimize -q 
%             subject to 
%             
%             for k = 1 : N1 - 1
%                 X( : , : , k) == hermitian_semidefinite(M_t) ; 
%                 diag(X( : , : , k)) >= 0
%                 for t = 1 : k 
%                     1 * q * noise_power +   1 * q * real(trace( H( : , : , k )*...
%                         sum( X(: , : , t + 1 : N1) , 3 )  ) - ...
%                             1 * trace(H( : , : , k ) * X(: , : , t))) <= 0
%                 end
%             end
%             k = N1 ;
%             X( : , : , k) == hermitian_semidefinite(M_t) ; 
%             diag(X( : , : , k)) >= 0
%             for t = 1 : k - 1 
%                 1 * q * noise_power +...
%                 1 * q * real(trace(  H( : , : , k )* sum( X(: , : , t + 1 : N1) , 3 )  )) - ...
%                 1 * real(trace(H( : , : , k ) * X(: , : , t))) <= 0
%             end
%             1 * q * noise_power - 1 * real(trace(H( : , : , k ) * X(: , : , k))) <= 0
%             trace(sum(X, 3)) <= P_T
%             
%         cvx_end
%         
%         if cvx_status( 1 : 6 ) == 'Solved'
% %             disp('It is feasible')
%             X_star = X ;
%             q_min = q ;
%         else
%             q_max = q ;
%         end
%     end
%     for k = 1 : N1
%             W( : , k) = GR(X_star( : , :,k) , 2) ;
%     end
%     q = q_min ;
%     rate = log2(q) ;
% X = X_star ;
% 
% % Evaluation of example 2 without GR
% r1 = log2(1 + (1 * trace(H( : , : , 1 ) * X(: , : , 1))) /...
%     (noise_power+ real(trace( H( : , : , 1 )*...
%     sum( X(: , : , 2) , 3 )  ) ))) ;
% r21 = log2(1 + (1 * trace(H( : , : , 2 ) * X(: , : , 1))) /...
%     (noise_power+ real(trace( H( : , : , 2 )*...
%     sum( X(: , : , 2) , 3 )  ) ))) ;
% r22 = log2(1 + (1 * trace(H( : , : , 2 ) * X(: , : , 2))) / noise_power) ;
% disp(r1)
% disp(r21)
% disp(r22)
% 
% % Evaluation of example 2 after GR
% r1 = log2(1 + (abs(h(:, 1)' * W(: , 1))^2)/(noise_power + abs(h(:, 1)' * W(: , 2))^2 )) ;
% r21 = log2(1 + (abs(h(:, 2)' * W(: , 1))^2)/(noise_power + abs(h(:, 2)' * W(: , 2))^2 )) ; 
% r22 = log2(1 + (abs(h(:, 2)' * W(: , 2))^2)/noise_power) ; 
% disp(r1)
% disp(r21)
% disp(r22)



