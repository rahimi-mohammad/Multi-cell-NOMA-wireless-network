% ------------------------------------------------------------------------
% rahimi-mohammad - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% PotentialGame - Computes the game table containing the utility function
% of the BS for each strategy.

% This function simulates a potential game for a scenario involving two
% base stations (BS1 and BS2) and multiple users. It computes the average
% achievable rate for each BS's user for different strategies of allocating
% IRS elements.

% Inputs:
%   M_t          - Number of transmitter antennas,
%   M_r          - Number of receiver antennas,
%   N1           - Number of users of BS1,
%   N2           - Number of users of BS2,
%   N_i          - Number of IRS elements,
%   P_T          - BS power,
%   x,y,z        - Users location,
%   x1,y1,z1     - BS1 location,
%   x2,y2,z2     - BS2 location,
%   x_i,y_i,z_i  - IRS location,
%   alpha_d      - BS to user pathloss,
%   alpha_r      - BS to IRS pathloss,
%   noise_power  - Noise power,
%   N_iter       - Number of iterations,
%   R            - User area center radius,
%   r            - Quantization step-size,
%   q_t          - Quantization step-size,
%   interference - Interference caused by IRS elements allocated to the
%                  other BS,
%   epsilon      - Parameter for rate calculation.

% Outputs:
%   utility - Utility function table,
%   rate    - Average achievable rate for each BS's user.

% ------------------------------------------------------------------------

function [utility, rate] = PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
                            x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
                            R, r, q_t, interference, epsilon)
    % Initialize utility and rate matrices
    rate = zeros(ceil(1/q_t) + 1, ceil(1/q_t) + 1, 2);
    utility = zeros(ceil(1/q_t) + 1, ceil(1/q_t) + 1, 2);
    
    % Initial user locations
    x0 = x;
    y0 = y;
    
    % Simulation loop
    for j = 1:N_iter
        % Generate random user locations
        t = 2 * pi * rand(N1, 1);
        radius = R(1) * sqrt(rand(N1, 1));
        x(1:N1) = x0(1:N1) + radius .* cos(t);
        y(1:N1) = y0(1:N1) + radius .* sin(t);
        
        t = 2 * pi * rand(N2, 1);
        radius = R(2) * sqrt(rand(N2, 1));
        x(N1 + 1:N1 + N2) = x0(N1 + 1:N1 + N2) + radius .* cos(t);
        y(N1 + 1:N1 + N2) = y0(N1 + 1:N1 + N2) + radius .* sin(t);
        
        % Compute channel gains
        [h_d1, G1, h_r1] = ChannelGain(M_t, M_r, N1, N_i, x(1:N1), y(1:N1), z(1:N1), x1, y1,...
                                                z1, x_i, y_i, z_i, alpha_d, alpha_r);
        [h_d2, G2, h_r2] = ChannelGain(M_t, M_r, N2, N_i, x(N1 + 1:N1 + N2), y(N1 + 1:N1 + N2), z(N1 + 1:N1 + N2),...
                                    x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r);
        
        % Perform game simulation based on the parameters
        if M_t == 1
            % Single antenna base station
            for s1 = 1:ceil(1/q_t) + 1
                for s2 = 1:ceil(1/q_t) + 2 - s1
                    % Number of IRS elements allocated to each BS
                    N_i1 = floor((s1 - 1) * N_i * q_t);
                    N_i2 = floor((s2 - 1) * N_i * q_t);
                    
                    % Compute rates for each user
                    [exact_rate1, ~] = Rate(h_d1, G1(1:N_i1), h_r1(1:N_i1, :), N1, N_i1, P_T, noise_power, 0, false, "NOMA");
                    [exact_rate2, ~] = Rate(h_d2, G2(1:N_i2), h_r2(1:N_i2, :), N2, N_i2, P_T, noise_power, 0, false, "NOMA");
                    
                    % Update cumulative rates
                    rate(s1, s2, 1) = rate(s1, s2, 1) + exact_rate1;
                    rate(s1, s2, 2) = rate(s1, s2, 2) + exact_rate2;
                end
            end
        elseif M_t > 1
            % Multiple antennas at the base station
            if interference == 0
                rate1 = zeros(ceil(1/q_t) + 1, 1);
                rate2 = zeros(ceil(1/q_t) + 1, 1);
                
                % Calculate rates for each strategy
                for s1 = 1:ceil(1/q_t) + 1
                    N_i1 = floor((s1 - 1) * N_i * q_t); % No. of elements allocated to BS1
                    [~, rate1(s1), ~] = Rate3(h_d1, G1(1:N_i1, :), h_r1(1:N_i1, :), N1, N_i1, P_T, noise_power, epsilon);
                end
                for s2 = 2:ceil(1/q_t) + 1
                    N_i2 = floor((s2 - 1) * N_i * q_t); % No. of elements allocated to BS2
                    [~, rate2(s2), ~] = Rate3(h_d2, G2(1:N_i2, :), h_r2(1:N_i2, :), N2, N_i2, P_T, noise_power, epsilon);
                end
                
                % Update cumulative rates
                for s1 = 1:ceil(1/q_t) + 1
                    for s2 = 1:ceil(1/q_t) + 2 - s1
                        rate(s1, s2, 1) = rate(s1, s2, 1) + rate1(s1);
                        rate(s1, s2, 2) = rate(s1, s2, 2) + rate2(s2);
                    end
                end
            end
        end
    end
    
    % Compute average rates
    rate = rate / N_iter;
end
