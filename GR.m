function y=GR(Z)
    N_r=100;
    N=size(Z,1)-1;
    [U,S,V]=svd(conj(Z));
%     z_hat=sqrt(S(1,1)*U(:,1,1)) ;
    V_prime=V*sqrt(S);
    z_hat=zeros(N+1,1);
    for i=1:N_r
        r=sqrt(0.5)*(randn(N+1,1)+1j*randn(N+1,1)); % r~CSCG(0,I)
        z_hat=z_hat+V_prime*r;
    end
    z_hat=z_hat/N_r;
    z_hat=z_hat./z_hat(N+1);
    abs(z_hat);
    angle(z_hat(1:N) );
    z=exp( 1j*angle(z_hat(1:N) ));
    y=z;
    
end