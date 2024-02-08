% ------------------------------------------------------------------------ 
% rahimi-mohammad - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% central_opt.m - This method returns achievable rates of the users of each 
% BS by solving a centralized optimization problem.
%-------------------------------------------------------------------------
% Inputs:
%   M_t                - No. transmitter antennas
%   M_r                - No. receiver antennas
%   N1                 - No. users of BS1,
%   N2                 - No. users of BS2,
%   N_i                - No. IRS elements,
%   P_T                - BS power,
%   x,y,z              - Users location,
%   x1,y1,z1           - BS1 location,   
%   x2,y2,z2           - BS2 location,   
%   x_i,y_i,z_i        - IRS location,      
%   alpha_d            - BS to user pathloss
%   alpha_r            - BS to IRS pathloss
%   noise_power        - Noise power
%   N_iter             - No. iterations
%   R                  - User area center radius
%   discrete           - Flag for discrete RIS phase shift
%   interference       - Flag for interference caused by IRS elements allocated to the other BS
%   epsilon            - Convergence criterion for optimization
%
% Outputs:
%   rate               - Average achievable rate for each BS's user (2x1).
%   Z_cvx              - IRS optimal phase shift ([e^(j*theta1) ... e^(j*thetaN)]).
% ------------------------------------------------------------------------

function [rate, Z_cvx] = central_opt(M_t, M_r, N1, N2, N_i, P_T, x, y, z,...
                                    x1, y1, z1, x2, y2, z2, x_i, y_i, z_i,...
                                    alpha_d, alpha_r, noise_power, N_iter,...
                                    R, discrete, interference, epsilon)
    rate = zeros(2, 1);
    x0 = x;
    y0 = y;
    N = N_i;

    if interference == 1
        error("RIS configuration sub-problem not implemented.");
        % Code for interference case
    else
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

            % Calculate channel gains
            [h_d1, G1, h_r1] = ChannelGain(M_t, M_r, N1, N_i, x(1:N1), y(1:N1), z(1:N1),...
                                            x1, y1, z1, x_i, y_i, z_i, alpha_d, alpha_r);
            [h_d2, G2, h_r2] = ChannelGain(M_t, M_r, N2, N_i, x(N1 + 1:N1 + N2), y(N1 + 1:N1 + N2),...
                                            z(N1 + 1:N1 + N2), x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r);

            if M_t == 1
                % RIS configuration sub-problem
                % Problem coefficients
                phi1 = zeros(N, N1); 
                Q1 = zeros(N + 1, N + 1, N1); 
                for i = 1:N1 
                    phi1(:, i) = diag(h_r1(:, i)') * G1; 
                    Q1(:, :, i) = [phi1(:, i) * phi1(:, i)' conj(h_d1(i)) * phi1(:, i);...
                                    h_d1(i) * phi1(:, i)' 0]; 
                end
                phi2 = zeros(N, N2); 
                Q2 = zeros(N + 1, N + 1, N2); 
                for i = 1:N2 
                    phi2(:, i) = diag(h_r2(:, i)') * G2; 
                    Q2(:, :, i) = [phi2(:, i) * phi2(:, i)' conj(h_d2(i)) * phi2(:, i);...
                                    h_d2(i) * phi2(:, i)' 0]; 
                end

                % Problem solver
                Q1 = 10^13 * Q1;
                Q2 = 10^13 * Q2;
                cvx_begin sdp quiet
                    variable X(N + 1, N + 1) complex hermitian; 
                    variable s(1, 1) complex hermitian; 
                    maximize s
                    subject to 
                    for m = 1:N1
                        10^13 * abs(h_d1(m))^2 + real(trace(Q1(:, :, m) * X)) >= s
                    end
                    for m = 1:N2
                        10^13 * abs(h_d2(m))^2 + real(trace(Q2(:, :, m) * X)) >= s
                    end

                    diag(X) == 1; 
                    X == hermitian_semidefinite(N + 1); 
                cvx_end

                Z_cvx = GR(X, 1);
                if discrete == true
                    % Apply discrete RIS phase shift
                    Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, 1);
                end

                Theta = diag(Z_cvx);

                % Power allocation sub-problem
                channel_gain1 = sort((abs(h_d1 + h_r1' * Theta * G1)).^2 * P_T / noise_power); 
                rate(1) = log2(1 + channel_gain1(end)) / N1;

                channel_gain2 = sort((abs(h_d2 + h_r2' * Theta * G2)).^2 * P_T / noise_power); 
                rate(2) = log2(1 + channel_gain2(end)) / N2;
            else
                error("Not implemented for multi-antenna basestation");
            end
        end
    end
end
