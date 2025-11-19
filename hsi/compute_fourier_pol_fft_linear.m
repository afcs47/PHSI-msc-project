function [S0, S1, S2, DoLP_fft, AoLP_fft] = compute_fourier_pol_fft_linear(I_theta, numAngles)% Compute Stokes parameters (linear only) using FFT

    % Reshape into [pixels x numAngles]
    [rows, cols, ~] = size(I_theta);
    I_pixels = reshape(I_theta, [], numAngles);

    % Apply FFT along angle dimension
    F = fft(I_pixels, [], 2) / numAngles; % normalized

    % Extract coefficients:
    % DC term - S0
    % cos(2θ), sin(2θ) - S1, S2
    a0 = 2 * real(F(:,1));   % DC component
    c2 = 2 * real(F(:,2));   % cos(2θ)
    s2 = -2 * imag(F(:,2));  % sin(2θ)

    % Map to Stokes parameters
    S0 = a0;
    S1 = c2;
    S2 = s2;

    % Derived polarization parameters
    DoLP_fft = sqrt(S1.^2 + S2.^2) ./ S0;
    AoLP_fft = 0.5 * atan2(S2, S1);

    % Reshape back to image
    S0 = reshape(S0, rows, cols);
    S1 = reshape(S1, rows, cols);
    S2 = reshape(S2, rows, cols);
    DoLP_fft = reshape(DoLP_fft, rows, cols);
    AoLP_fft = reshape(AoLP_fft, rows, cols);
end
