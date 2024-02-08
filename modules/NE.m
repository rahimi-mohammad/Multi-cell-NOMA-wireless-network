% ------------------------------------------------------------------------
% rahimi-mohammad - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% NE - Calculates Nash Equilibrium for a two-player game with two pure
% strategies.

% This function computes the probability of each player choosing the "W/O"
% strategy to achieve Nash Equilibrium.

% Inputs:
%   u - Utility function of players. It's a 3D array representing the
%       utility of each player for each combination of strategies.

% Outputs:
%   P - Probability of playing strategy "W/O" for each player.
% ------------------------------------------------------------------------

function P = NE(u)
    % Initialize the probability vector
    P = zeros(2, 1);
    
    % Check for dominant strategies and compute probabilities accordingly
    if (u(1, 1, 1) - u(2, 1, 1)) * (u(1, 2, 1) - u(2, 2, 1)) < 0
        % Player 1 has no dominant strategy
        if (u(1, 1, 2) - u(1, 2, 2)) * (u(2, 1, 2) - u(2, 2, 2)) < 0
            % Player 2 has no dominant strategy
            P(1) = (u(2, 2, 2) - u(2, 1, 2)) / (u(1, 1, 2) + u(2, 2, 2) - u(1, 2, 2) - u(2, 1, 2));
            P(2) = (u(1, 2, 1) - u(2, 2, 1)) / (u(1, 2, 1) - u(2, 2, 1) + u(2, 1, 1) - u(1, 1, 1));
        elseif u(1, 1, 2) > u(1, 2, 2)
            P(2) = 1;
            P(1) = (1 + sign(u(1, 1, 1) - u(2, 1, 1))) / 2;
        else
            P(1) = (1 + sign(u(1, 2, 1) - u(2, 2, 1))) / 2;
        end
    elseif u(1, 1, 1) > u(2, 1, 1)
        P(1) = 1;
        P(2) = (1 + sign(u(1, 1, 2) - u(1, 2, 2))) / 2;
    else
        P(2) = (1 + sign(u(2, 1, 2) - u(2, 2, 2))) / 2;
    end
end
