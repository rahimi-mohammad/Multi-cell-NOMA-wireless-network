% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method calculates aceivable rate of the BS's users for each
% strategy
% Inputs:
    % h_d                - BS to users channel gain,
    % G                  - BS to IRS channel gain,
    % h_r                - IRS to users channel gain.    
    % N1                 - No.  users of BS1,
    % N                  - No.  IRS elements(It can be a vector),
    % noise_power        - Noise power
% Outputs:
    % rate               - Acheivable Rate of the BS(the lower bound).
    % Z_cvx              - IRS optimal Phase shift.
% ------------------------------------------------------------------------
function [W , rate , Z_cvx] = Rate3(h_d, G, h_r, N1, N, P_T, noise_power , epsilon)
    M_t = size(h_d, 1) ;
    W = zeros(M_t , N1) ;   
    W = 1 * randn(M_t , N1)/sqrt(4 * N1 * M_t/ P_T) ; 
%     load('W.mat')
    Z_cvx = zeros(N + 1, 1) ; 
    rate = 0 ; 
    if N == 0
        % without IRS
        Z_cvx = [] ; 
        if M_t==1
            h_d = sort(abs(h_d)) ; 
        rate = log2(1 + P_T * h_d(1)^2/noise_power)/N1 + rate ;  
        elseif M_t>1
            % I will write this later
        end 
    else
        % with IRS
        if M_t == 1
            %% problem coefficients
            phi = zeros(N, N1) ; 
            Q = zeros(N + 1, N + 1, N1) ; 
            for i = 1:N1 
                phi(:,i) = diag(h_r(:,i)') * G ; 
                Q(:,:,i) = [phi(:,i) * transpose(conj(phi(: , i))) conj(h_d(i)) * phi(:,i) ; h_d(i) * transpose(conj(phi(:,i))) 0] ; 
            end
            %% problem solver
            Q = 10^13 * Q ; 
            cvx_begin sdp quiet
                variable X(N + 1,N + 1) complex hermitian  ;
                variable s(1, 1)     complex hermitian  ;
                maximize s
                subject to 
                for m = 1 : N1
                    10^13 * abs(h_d(m))^2 + real(trace(Q(:, :, m) * X)) >= s
                end
                diag(X) == 1 ; 
                X == hermitian_semidefinite(N + 1) ; 
            cvx_end
            Z_cvx = GR(X) ;
        %                 Q=10^(-13) * Q ; 
            %% channel gains
            channel_gain=sort(abs(h_d + h_r' * diag(Z_cvx) * G)) ; 
            rate_lower_bound=log2(1 + P_T * channel_gain(1)^2/noise_power)/N1 ; 
            rate=rate + rate_lower_bound ; 
        elseif M_t>1
            % I will write this later
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve P1.2 IRS optimization
            [h_d, G, h_r] = CCSOrder(h_d, G, h_r) ;
            h = zeros(M_t , N1) ; % Cascade channel
            rate0 = 0 ;
            rate2 = 0.2 ;
            feasible = false ;
            while feasible == false
                disp('regenerating W') 
                W = 1 * randn(M_t , N1)/sqrt(4 * N1 * M_t/ P_T) ; 
                [rate1, Z_cvx , feasible] = MISONoma1(W , epsilon, h_d , G , h_r , noise_power) ;
            end
            while abs(rate2 - rate0) > 0.1
            rate0 = rate2 ;
                for k = 1 : N1
                    h( : , k ) = (h_r(: , k)' * diag(Z_cvx) * G + h_d(: , k)')' ;
                end
            %% Solve P1.4
%             Z_cvx
            [rate2, W] = MISONoma2(epsilon, h , noise_power , P_T) ;
%             disp(['step1 : ' , num2str(rate1) , ' , step2 : ' , num2str(rate2)])
            [rate1, Z_cvx , feasible] = MISONoma1(W , epsilon, h_d , G , h_r , noise_power) ;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            rate = rate2 ;

            %% problem coefficients
%             phi = zeros(N, N1) ; 
%             Q = zeros(N + 1, N + 1, N1) ; 
%             for i = 1:N1
%                 phi(:, i) = diag(h_r(:, i)') * G ; 
%                 Q(:, :, i) = [phi(:, i) * transpose(conj(phi(:, i))) conj(h_d(i)) * phi(:, i) ; h_d(i) * transpose(conj(phi(:, i))) 0] ; 
%             end
        end
    end

end