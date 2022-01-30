
cvx_begin sdp
    variable X(N+1,N+1) complex hermitian  ;
    
    maximize real(trace(Q(:,:,1)*X))    %% trace(A4*X) is Real if X and A4 is hemitian matrix
    subject to 
     diag(X) == 1;
     X==hermitian_semidefinite(N+1);
cvx_end