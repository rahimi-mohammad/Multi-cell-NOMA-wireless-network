% ------------------------------------------------------------------------
% rahimi-mohammad - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method calculates RIS discrete phase shift based on the
% norm.
% Inputs:
%   Z_cvx        - RIS continuous phase shift
%   N_bits       - No. of RIS phase shift discretization bits
%   method       - The method used for conversion (optional)
% Outputs:
%   Z_cvx        - RIS discrete phase shift
% ------------------------------------------------------------------------

function Z_cvx = discrete_RIS_phase_shift(Z_cvx, N_bits, method)
    % Calculate the discrete phase shift
    Z_discrete = exp(1j * 2 * pi * round((pi + angle(Z_cvx)) / 2 / pi * 2^N_bits) / 2^N_bits);

    % Update the continuous phase shift with the discrete values
    Z_cvx = Z_discrete;
end
