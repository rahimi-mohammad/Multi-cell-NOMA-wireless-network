function P=NE(u)
%This function calculates Nash Equilibrium for a two-player game with two pure strategies. 
    
    %               Game Table:
 %  ____________________________________________________________________
% |    |            with probability P(2)     |   1-P(2)               |
% |    |        W/O                           |     W                  |
% |___________________________________________|________________________|
% | W/O|    u(1,1,1),u(1,1,2)                 |    u(1,2,1),u(1,2,2)   |
% | W  |    u(2,1,1),u(2,1,2)                 |    u(2,2,1),u(2,2,2)   |
% |____|______________________________________|________________________|
% Inputs:
    % u               % Utility function of players.
% Outputs:
    % P                - Probability of playing strategy W/O.
%-------------------------------------------------------------------------
    
    P=zeros(2,1);
    if (u(1,1,1)-u(2,1,1))*(u(1,2,1)-u(2,2,1))<0         % user 1 has no dominant strategy
        if (u(1,1,2)-u(1,2,2))*(u(2,1,2)-u(2,2,2))<0     % user 2 has no dominant strategy 
            P(1)=(u(2,2,2)-u(2,1,2))/(u(1,1,2)+u(2,2,2)-u(1,2,2)-u(2,1,2));
            P(2)=(u(1,2,1)-u(2,2,1))/(u(1,2,1)-u(2,2,1)+u(2,1,1)-u(1,1,1));
        elseif u(1,1,2)>u(1,2,2)
            P(2)=1;
            P(1)=(1+sign(u(1,1,1)-u(2,1,1)))/2;
        else
            P(1)=(1+sign(u(1,2,1)-u(2,2,1)))/2;
        end
    elseif u(1,1,1)>u(2,1,1)
        P(1)=1;
        P(2)=(1+sign(u(1,1,2)-u(1,2,2)))/2;
    else
        P(2)=(1+sign(u(2,1,2)-u(2,2,2)))/2;
    end
end