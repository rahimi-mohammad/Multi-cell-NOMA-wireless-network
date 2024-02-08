% ------------------------------------------------------------------------
% rahimi-mohammad - Sharif University of Technology, Iran
% GR.m - Generate RIS phase shift using different methods.
% ------------------------------------------------------------------------
% This function generates RIS phase shift based on different methods:
% 1. Random method: Averages multiple random vectors in the direction of the
%    principal eigenvector of the input matrix.
% 2. Minimum norm method: Finds the vector with the minimum Frobenius norm
% difference among vectors that result in the same matrix as the input.
%
% Inputs:
%   Z       - Input matrix representing RIS channel gains
%   method  - Method for generating RIS phase shift (1 for Random, 2 for
%             Minimum Norm)
%
% Outputs:
%   y       - RIS phase shift vector generated using the specified method
% ------------------------------------------------------------------------

function y = GR(Z, method)
    N_r = 1000; % Number of random trials
    
    % Random method
    if method == 1  
        N = size(Z, 1) - 1; 
        [~, S, V] = svd(conj(Z)); 
        V_prime = V * sqrt(S); 
        z_hat = zeros(N + 1, 1); 
        for i = 1:N_r
            r = sqrt(0.5) * (randn(N + 1, 1) + 1j * randn(N + 1, 1));  % r~CSCG(0,I)
            z_hat = z_hat + V_prime * r; 
        end
        z_hat = z_hat / N_r; 
        z_hat = z_hat ./ z_hat(N + 1); 
        z = exp(1j * angle(z_hat(1 : N))); 
        y = z; 

    % Minimum norm method
    elseif method == 2
        N = size(Z, 1); 
        [~, S, V] = svd(conj(Z)); 
        V_prime = V * sqrt(S); 
        z_hat = zeros(N, 1); 
        for i = 1:N_r
            r = sqrt(0.5) * (randn(N, 1) + 1j * randn(N, 1));  % r~CSCG(0,I)
            if norm(z_hat * z_hat' - Z) > norm(V_prime * r * (V_prime * r)' - Z) 
                z_hat = V_prime * r;
            end
        end
        y = z_hat; 
    end
end
