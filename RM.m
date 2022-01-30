function y=RM(Z)
   epsilon=1e-2;
    X=Z;
    N=size(X,2);
    delta=ones(N,N);
    while delta~=0
        [U,S,V]=svd(X);
         V=V*sqrt(S);
         cvx_begin sdp
            variable delta(N,N) complex hermitian  ;
            maximize 1       
            subject to 
                for i=1:N
                   A=zeros(N,N);
                   A(i,i)=1;
                   trace(A*V*delta*V')==1; 
                end
            delta==hermitian_semidefinite(N);
         cvx_end
         [U,S,~]=svd(delta);
        alpha=-1/S(1,1)
        X=V*(eye(N)+alpha*delta)*V';
         disp('|delta|=');
         disp(num2str(abs(delta)));
    end
    y=V(1:end-1,1);
end