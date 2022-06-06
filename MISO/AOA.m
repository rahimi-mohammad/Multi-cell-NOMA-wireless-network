% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% ChannelGain.m - This method generates channel gains of the scenario
% Inputs:
    % h_d               % Direct channel gain,
    % G                 % IRS-BS channel gain
    % h_r               % IRS-users channel gain,  
    % N1                 - No.  users of BS1,
    % N_i                - No.  IRS elements,
    % P_T                - BS power,
    % x,y,z              - Users location,
    % x1,y1,z1           - BS location,   
    % x_i,y_i,z_i        - IRS location,      
    % alpha_d            - BS to user pathloss,
    % alpha_r            - BS to IRS pathloss.
    % noise_power            
% Outputs:
    % theta                - Scenario: s=0 without IRS and s=1 with RIS,
    % G                  - No.  users of BS1,
    % h_r                - No.  IRS elements.    
% ------------------------------------------------------------------------

function theta=AOA(h_d, G, h_r, w)
h_d=transpose(h_d)
N1=size(h_d,2)
N=size(G,1)
M_t=size(h_d,1)
h_dt=zeros(N1,1)
gt=zeros(N,1)

for i=1:N1
   h_dt(i)=h_d(:,i)'*w(:,i); 
   G(:,i)
   gt(i)=G(:,i)'*w(:,i);     

end
h_dt 
gt(i)
   

% 
%     while cconverge!=true
%        % Solve first subproblem 
% 
% 
% 
% 
% 
%        % Solve second subproblem 
%     end


end