% ------------------------------------------------------------------------
% rahimi-mohammad - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate - Calculates achievable rate of the BS's users for each strategy.
%
% This function computes the achievable rate of the BS's users for each
% strategy considering different MAC protocols and the presence of an IRS.
%
% Inputs:
%   h_d          - BS to users channel gain,
%   G            - BS to IRS channel gain,
%   h_r          - IRS to users channel gain,
%   N1           - Number of users of BS1,
%   N_i          - Number of IRS elements (It can be a vector),
%   P_T          - BS power,
%   noise_power  - Noise power,
%   N_bits       - No. RIS phase shift discretization bits,
%   discrete     - Continuous or discrete RIS phase shift,
%   MAC          - Multiple Access Control.
%
% Outputs:
%   rate         - Average achievable rate of the BS (EACH USER).
%   Z_cvx        - IRS optimal Phase shift ([e^(j*theta1) ... e^(j*thetaN)]).
%
% ------------------------------------------------------------------------

function [rate, Z_cvx] = Rate(h_d, G, h_r, N1, N_i, P_T, noise_power, N_bits, discrete, MAC)
    M_t = size(h_d, 2);

    if N_i == 0
        % Without IRS
        Z_cvx = [];

        if M_t == 1 % SISO
            h_d = sort(abs(h_d));
            switch MAC
                case 'NOMA'
                    if N1 == 2
                        % NOMA specific calculation
                        channel_gain = h_d.^2 * P_T / noise_power;
                        q_star = (channel_gain(1) - channel_gain(2) + ...
                            sqrt((channel_gain(1) - channel_gain(2))^2 + ...
                            4 * channel_gain(1) * (channel_gain(2) + ...
                            channel_gain(1) * channel_gain(2)))) / 2 / channel_gain(1);
                        rate = log2(q_star);
                    else
                        % General NOMA calculation
                        [rate, ~] = SISONoma(0.1, h_d.^2 * P_T / noise_power);
                    end
                case 'FDMA'
                    % FDMA calculation
                    q_star = (max(abs(h_d))^2 * P_T / noise_power) / N1;
                    rate = log2(1 + q_star);
                case 'TDMA'
                    % TDMA calculation
                    TDMA_rate = log2(1 + P_T / noise_power * h_d.^2) / N1;
                    rate = sum(TDMA_rate) / N1;
                otherwise
                    error('Invalid MAC');
            end
        elseif M_t > 1
            error('Not implemented');
        end
    else
        % With IRS
        [rate, Z_cvx] = RIS_configuration(h_d, G, h_r, N1, N_i, P_T, ...
            noise_power, M_t, N_bits, discrete, MAC);
    end
end
