function y = GR(Z , method)
    N_r = 500 ; 
if method == 1
    N = size(Z,1) - 1 ; 
    [~,S,V] = svd(conj(Z)) ; 
    V_prime = V * sqrt(S) ; 
    z_hat = zeros(N + 1,1) ; 
    for i = 1 : N_r
        r = sqrt(0.5) * (randn(N + 1,1) + 1j * randn(N + 1,1)) ;  % r~CSCG(0,I)
        z_hat = z_hat + V_prime * r ; 
    end
    z_hat = z_hat/N_r ; 
    z_hat = z_hat./z_hat(N + 1) ; 
    abs(z_hat) ; 
    angle(z_hat(1 : N) ) ; 
    z = exp( 1j * angle(z_hat(1 : N) )) ; 
    y = z ; 
    
    
    
elseif method == 2
        N = size(Z,1) ; 
    [~,S,V] = svd(conj(Z)) ; 
    V_prime = V * sqrt(S) ; 
    z_hat = zeros(N,1) ; 
    for i = 1 : N_r
        r = sqrt(0.5) * (randn(N ,1) + 1j * randn(N ,1)) ;  % r~CSCG(0,I)
        if norm(z_hat * z_hat' - Z) > norm(V_prime * r * (V_prime * r)' -Z) 
            z_hat = V_prime * r ;
        end
    end
 y = z_hat ; 

% elseif method == 3
%     mu = [0 0];
%     Sigma = Z ;
%     % rng('default')  % For reproducibility
%     R = mvnrnd(mu,Sigma,N_r);
%     for k = 1 : N_r
% %         R( k , :)/()
% 
%     end
end

end
