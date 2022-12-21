h_d=h_d1;
G=G1;
h_r=h_r1;
h = zeros(M_t , N1) ;

W = ones(M_t , N1)/N1 ; 
rate0 = 0 ;
rate = 0.2 ;

while abs(rate - rate0) > 0.1
    rate0 = rate 
    [rate, Z_cvx] = MISONoma1(W , epsilon, h_d , G , h_r , noise_power) 
    for k = 1 : N1
        h( : , k ) = (h_r(: , k)' * diag(Z_cvx) * G + h_d(: , k)')' ;
    end
    %% Solve P1.4
    [rate, W] = MISONoma2( epsilon, h , noise_power , P_T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end