% ------------------------------------------------------------------------
% mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% RIS_configuration - Calculates the RIS phase shift and achievable rate when
%                      an RIS is present in the setup.
%
% This function computes the RIS phase shift and achievable rate considering
% different MAC protocols and the presence of an RIS.
%
% Inputs:
%   h_d          - BS to users channel gain,
%   G            - BS to IRS channel gain,
%   h_r          - IRS to users channel gain.
%   N1           - No. users of BS1,
%   N_i          - No. IRS elements (It can be a vector),
%   P_T          - BS power,
%   noise_power  - Noise power,
%   M_t          - No. transmitter antennas,
%   N_bits       - No. RIS phase shift discretization bits,
%   discrete     - Continuous or discrete RIS phase shift,
%   MAC          - Multiple Access Control.
%
% Outputs:
%   rate         - Achievable Rate of the BS.
%   Z_cvx        - IRS optimal Phase shift([e^(j*theta1) ... e^(j*thetaN)]).
%
% ------------------------------------------------------------------------

function [rate, Z_cvx] = RIS_configuration(h_d, G, h_r, N1, N_i, P_T, noise_power, M_t, N_bits, discrete, MAC)
    N = N_i;
    rate = 0;

    switch MAC
        case 'NOMA'
            if M_t == 1 % Single-antenna
                if N1 == 2
                Z_cvx =  exp(1j*(angle(h_r(: , 2)) + angle(h_d(2)) - angle(G(:)))) ;
                    if discrete == true
                        Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, 1);
                    end
                    Theta = diag(Z_cvx);
                    channel_gain = sort(abs(h_d + h_r' * Theta * G)).^2 * P_T / noise_power;
                    q_star = (channel_gain(1) - channel_gain(2) + ...
                        sqrt((channel_gain(1) - channel_gain(2))^2 + ...
                        4 * channel_gain(1) * (channel_gain(2) + ...
                        channel_gain(1) * channel_gain(2)))) / 2 / channel_gain(1);
                    rate = log2(q_star);
                else
                    % Problem coefficients
                    phi = zeros(N, N1);
                    Q = zeros(N + 1, N + 1, N1);
                    for i = 1:N1
                        phi(:, i) = diag(h_r(:, i)') * G;
                        Q(:, :, i) = [phi(:, i) * transpose(conj(phi(:, i))) conj(h_d(i)) * phi(:, i); ...
                                      h_d(i) * transpose(conj(phi(:, i))) 0];
                    end
                    % Problem solver
                    Q = 10^13 * Q;
                    cvx_begin sdp quiet
                        variable X(N + 1, N + 1) complex hermitian;
                        variable s(1, 1) complex hermitian;
                        maximize s
                        subject to
                            for m = 1:N1
                                10^13 * abs(h_d(m))^2 + real(trace(Q(:, :, m) * X)) >= s
                            end
                            diag(X) == 1;
                            X == hermitian_semidefinite(N + 1);
                    cvx_end
                    Z_cvx = GR(X, 1);
                    if discrete == true
                        Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, 1);
                    end
                    Q = 10^(-13) * Q;
                    Theta = diag(Z_cvx);
                    channel_gain = sort(abs(h_d + h_r' * Theta * G));
                    rate = log2(1 + P_T * channel_gain(end)^2 / noise_power) / N1;
                end
            elseif M_t > 1
                % Not implemented yet
                error('Use Rate3 for MISO scenario');
            end

        case 'FDMA'
            if M_t == 1 % Single-antenna
                if N1 == 2
                    Z_cvx = exp(1j * (angle(h_r(:, 2)) + angle(h_d(2)) - angle(G(:))));

                    if discrete == true
                        Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, 1);
                    end
                    Theta = diag(Z_cvx);
                    channel_gain = sort(abs(h_d + h_r' * Theta * G)).^2 * P_T / noise_power;
                    q_star = (channel_gain(1) - channel_gain(2) + ...
                        sqrt((channel_gain(1) - channel_gain(2))^2 + ...
                        4 * channel_gain(1) * (channel_gain(2) + ...
                        channel_gain(1) * channel_gain(2)))) / 2 / channel_gain(1);
                    rate = log2(q_star);
                else
                    % Problem coefficients
                    phi = zeros(N, N1);
                    Q = zeros(N + 1, N + 1, N1);
                    for i = 1:N1
                        phi(:, i) = diag(h_r(:, i)') * G;
                        Q(:, :, i) = [phi(:, i) * transpose(conj(phi(:, i))) conj(h_d(i)) * phi(:, i); ...
                                      h_d(i) * transpose(conj(phi(:, i))) 0];
                    end
                    % Problem solver
                    Q = 10^13 * Q;
                    cvx_begin sdp quiet
                        variable X(N + 1, N + 1) complex hermitian;
                        variable s(1, 1) complex hermitian;
                        maximize s
                        subject to
                            for m = 1:N1
                                10^13 * abs(h_d(m))^2 + real(trace(Q(:, :, m) * X)) >= s
                            end
                            diag(X) == 1;
                            X == hermitian_semidefinite(N + 1);
                    cvx_end
                    Z_cvx = GR(X, 1);
                    if discrete == true
                        Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, 1);
                    end
                    Q = 10^(-13) * Q;
                    Theta = diag(Z_cvx);
                    channel_gain = sort(abs(h_d + h_r' * Theta * G));
                    q_star = (max(abs(channel_gain))^2 * P_T / noise_power) / N1;
                    rate = log2(1 + q_star) / N1; % FDMA
                end
            elseif M_t > 1
                % Not implemented yet
                error('Use Rate3 for MISO scenario');
            end

        case 'TDMA'
            TDMA_rate = 1 / N1 * log2(1 + P_T / noise_power * (abs(h_d) + abs(h_r') * abs(G)).^2);
            rate = sum(TDMA_rate) / N1;
            Z_cvx = exp(1j * (angle(h_d) - angle(h_r') + angle(G')));
            if discrete == true
                Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, 1);
            end

        otherwise
            error('Invalid MAC');
    end
end
