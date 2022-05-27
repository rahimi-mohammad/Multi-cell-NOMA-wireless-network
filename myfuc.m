% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% M.m - This method calculates aceivable rate of the BS's users for each
% strategy
% Inputs:
%     s                  - Scenario: s=0 without IRS and s=1 with RIS,
    % N1                 - No.  users of BS1,
    % N_i                - No.  IRS elements,
    % P_T                - BS power,
    % x,y,z              - Users location,
    % x1,y1,z1           - BS location,   
    % x_i,y_i,z_i        - IRS location,      
    % alpha_d            - BS to user pathloss
    % alpha_r            - BS to IRS pathloss
    % noise_power        - Noise power
    % N_iter             - No. iteraions
% Outputs:
    % U                  - Acheivable Rate of the BS.
% ------------------------------------------------------------------------
function X=myfuc(epsilon)
    wd
    wd
    b0=;
    b=b0;
    a=0+epsilon;
    
    while b-a>epsilon
       
        if feasible=='True'
           a=b; 
           b=(b+b0)/2;
        else
           b=(a+b)/2; 
        end
            
        
        
        
    end
    Q=10^13*Q;
        cvx_begin sdp
            variable X(N+1,N+1) complex hermitian  ;
            variable s(1,1)     complex hermitian  ;
            maximize s
%             min(10^13*abs(h_d(1))^2+real(trace(Q(:,:,1)*X)),10^1e3*abs(h_d(2))^2+real(trace(Q(:,:,2)*X)))      %% trace(A4*X) is Real if X and A4 is hemitian matrix
            subject to 
            for m=1:N1
                10^13*abs(h_d(m))^2+real(trace(Q(:,:,m)*X))>=s
            end
            diag(X) == 1;
            X==hermitian_semidefinite(N+1);
        cvx_end

    
    
    
    
end