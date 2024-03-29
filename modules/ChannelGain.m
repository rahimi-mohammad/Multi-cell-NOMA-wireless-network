% ------------------------------------------------------------------------
% rahimi-mohammad - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% ChannelGain.m - This method generates channel gains of the scenario.
% The cascade channel for all users will be h_d + h_r^H * Theta * G
% Inputs:
%   M_t             - No. of transmitter antennas
%   M_r             - No. of receiver antennas
%   N1              - No. users of BS1
%   N               - No. IRS elements
%   x, y, z         - Users location
%   x1, y1, z1      - BS location
%   x_i, y_i, z_i   - IRS location
%   alpha_d         - BS to user pathloss
%   alpha_r         - BS to IRS pathloss
%   noise_power     - Noise power
% Outputs:
%   h_d             - BS to users channel gain (N1 * M_t)
%   G               - BS to IRS channel gain (N * M_t)
%   h_r             - IRS to users channel gain (N * N1)
% ------------------------------------------------------------------------

function [h_d, G, h_r] = ChannelGain(M_t, M_r, N1, N, x, y, z, x1, y1, z1,...
                                    x_i, y_i, z_i, alpha_d, alpha_r)
    % Generate BS to IRS channel gain
    G = (10^0.2) * (randn(N, M_t) + 1j * randn(N, M_t)) * sqrt(1/2) * ...
        ((10^0.2) / (sqrt((x_i - x1).^2 + (y_i - y1).^2 + (z_i - z1).^2))^alpha_r);
    
    % Generate BS to users channel gain
    h_d = (10^0.2) * (randn(N1, M_t) + 1j * randn(N1, M_t)) * sqrt(1/2) .* ...
        ((10^0.2) ./ sqrt((x - x1).^2 + (y - y1).^2 + (z - z1).^2).^alpha_d);
    
    % Generate IRS to users channel gain
    h_r = zeros(N, N1);
    for i = 1:N1
        h_r(:, i) = (10^0.2) * (randn(N, 1) + 1j * randn(N, 1)) * sqrt(1/2) .* ...
            ((10^0.2) / sqrt((x(i) - x_i).^2 + (y(i) - y_i).^2 + (z(i) - z_i).^2).^alpha_r);
    end
end
