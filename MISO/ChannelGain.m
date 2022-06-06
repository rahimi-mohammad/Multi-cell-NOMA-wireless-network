% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% ChannelGain.m - This method generates channel gains of the scenario
% Inputs:
    % M_t               % No. of transmitter antennas,
    % M_r               % No. of receiver antennas,   
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
    % h_d                - Scenario: s=0 without IRS and s=1 with RIS,
    % G                  - No.  users of BS1,
    % h_r                - No.  IRS elements.    
% ------------------------------------------------------------------------
function [h_d, G, h_r]=ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                    x_i, y_i, z_i, alpha_d, alpha_r)
    if M_r==1
        G=(10^0.2)*(randn(N,M_t,M_r)+1j*randn(N,M_t,M_r))*sqrt(1/2).*...
          ((10^0.2)./(sqrt( (x_i-x1).^2+(y_i-y1).^2+(z_i-z1).^2 )).^alpha_r);
        h_d=(10^0.2)*(randn(N1,M_t,M_r)+1j*randn(N1,M_t,M_r))*sqrt(1/2).*...
            ((10^0.2)./(sqrt( (x-x1).^2+(y-y1).^2+(z-z1).^2 )).^alpha_d);
        
        h_r=zeros(N,N1);
%         x
%         
%        y
%        z
%        
        for i=1:N1
            h_r(:,i)=(10^0.2)*(randn(N,M_r)+1j*randn(N,M_r))*sqrt(1/2).*...
                ((10^0.2)./(sqrt( (x(i)-x_i).^2+(y(i)-y_i).^2+(z(i)-z_i).^2 )).^alpha_r);
        end
    end
end