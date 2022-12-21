function [h_d, G, h_r] = CCSOrder(h_d, G, h_r) 
    N1 = size(h_d , 2) ; % No. users
    N = size(h_r , 1) ; % No. IRS elements
    M_t = size(h_d , 1) ; % No. transmitter antennas
    S = zeros( N + 1 , N + 1, N1 );
    J = zeros( N , M_t , N1) ;  %  N * M_t
    h = zeros(N1 , 1) ;
    h_d1 = zeros(size(h_d)) ;
    h_r1 = zeros(size(h_r)) ;
    for k = 1 : N1
        J( : , : , k ) = diag(h_r(: , k)') * G ; % N * M_t
        S( : , : ,k) =1e10 * [J(: , : , k) * J(: , : , k)' , J(: , : , k) * h_d(  : , k) ;...
            h_d( : , k )' * J(: , : , k)' , 0] ; 
    end
    for k = 1 : N1
%% problem solver
        cvx_begin sdp quiet 
            variable X(N + 1,N + 1) complex hermitian  ;
            minimize -1 *real(trace(S( : , : ,k) * X)) -...
                    1e10 *real(h_d(: , k)' * h_d(: , k))
            subject to 
            X == hermitian_semidefinite(N + 1) ; 
            diag(X) == 1 ; 
        cvx_end
        Z_cvx = GR(X , 1) ;
        h(k) = norm(h_d(: , k)' + h_r(: , k)' * diag(Z_cvx) * G ) ;
    end
    [~ , order] = sort(h , 'ascend') ;
    for k = 1 : N1
        h_d1(: ,k ) = h_d(: ,order(k) ) ;
        h_r1(: ,k ) = h_r(: , order(k) ) ;
    end
    h_d = h_d1 ;
    h_r = h_r1 ;
    