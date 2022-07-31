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
    % N_i                - No.  IRS elements(It can be a vector),
    % noise_power        - Noise power
    % N_iter             - No. iteraions
% Outputs:
    % rate               - Acheivable Rate of the BS(the lower bound).
    % Z_cvx              - IRS optimal Phase shift.
% ------------------------------------------------------------------------
function [rate, Z_cvx] = Rate2(h_d, G, h_r, N1, N_i, P_T, noise_power, N_iter)
    Z_cvx = zeros(N_i + 1, 1) ; 
    final_rate = 0 ; 
    N = N_i ; 
    M_t = size(h_d, 2) ; 
    if N_i==0
        % without IRS
        Z_cvx = [] ; 
        if M_t==1
            h_d = sort(abs(h_d)) ; 
        final_rate = log2(1 + P_T * h_d(1)^2/noise_power)/N1 + final_rate ;  
        elseif M_t>1
            % I will write this later
        end 
        final_rate = final_rate/N_iter ; 
    else
        % with IRS
        if M_t==1
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
            final_rate=final_rate + rate_lower_bound ; 
            final_rate=final_rate/N_iter ; 
        elseif M_t>1
            % I will write this later
            %% problem coefficients
            phi = zeros(N, N1) ; 
            Q = zeros(N + 1, N + 1, N1) ; 
            for i = 1:N1
                phi(:, i) = diag(h_r(:, i)') * G ; 
                Q(:, :, i) = [phi(:, i) * transpose(conj(phi(:, i))) conj(h_d(i)) * phi(:, i) ; h_d(i) * transpose(conj(phi(:, i))) 0] ; 
            end
        end
    end
    rate = final_rate ; 

end