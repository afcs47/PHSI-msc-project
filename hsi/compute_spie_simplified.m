function [S0, S1, S2, DoLP_img, AoP_img] = compute_spie_simplified(I_theta)

    [rows, cols, numAngles] = size(I_theta);
    I_reshaped = reshape(I_theta, [], numAngles); % [pixels, angles]
    
    % Perform FFT along angles
    F = fft(I_reshaped, [], 2) / numAngles;
    
    % Extract relevant components
    F0 = real(F(:, 1)); % DC component - S0 * C1/2
    F2 = F(:, 3); % Harmonic - S1, S2
    
    % Compute polarization parameters
    DoLP = abs(F2) ./ F0;                          % [pixels x 1]
    AoP_rad = 0.5 * atan2(imag(F2), real(F2));     % [pixels x 1]
    AoP_deg = rad2deg(AoP_rad);

    % Reshape to image
    DoLP_img = reshape(DoLP, rows, cols);
    AoP_img = reshape(AoP_deg, rows, cols);

    S0 = reshape(F0, rows, cols);
    S1 = reshape(real(F2), rows, cols);
    S2 = reshape(imag(F2), rows, cols);

end
