function [DoLP, AoLP] = calculate_polarization_fourier_4th(polar_images) % outputs DoLP, AoLP, computed from 4th order Fourier model from the inputed polar_images demosaiced polarization stack
    % Define the known polarization angles (in degrees) used for capturing the images
    angles_deg = [0, 45, 90, 135];

    % Create the matrix for 4th-order fit
    A4 = [ones(4,1), cosd(2 * angles_deg)', sind(2 * angles_deg)', cosd(4 * angles_deg)', sind(4 * angles_deg)']; % 4th-order model: [a0, a2, b2, a4, b4]

    % Initialize result arrays
    [H, W, C, ~] = size(polar_images);
    DoLP = zeros(H, W, C);
    AoLP = zeros(H, W, C);

    % Loop through each color channel separately
    for c = 1:C
        % Extract current color channel across all angles
        I = squeeze(polar_images(:,:,c,:)); % [H, W, 4]
        % Reshape for matrix operations: each pixel is a row -> shape [H*W, 4], transpose to [4, N] and convert to double
        I_reshaped = double(reshape(I, [], 4)'); % [4, H*W]

        % Solve linear least squares for 4th-order coefficients
        coeffs = A4 \ I_reshaped; % [5, H*W]
        a0 = coeffs(1, :);
        a2 = coeffs(2, :);
        b2 = coeffs(3, :);
        % In case a4, b4 may be used for advanced metrics, uncomment following line
        % a4 = coeffs(4, :); b4 = coeffs(5, :);

        % Compute DoLP and AoLP from Fourier coefficients
        dolp = sqrt(a2.^2 + b2.^2) ./ (a0 + eps);
        aolp = 0.5 * atan2(b2, a2);

        % Reshape back to [H, W]
        DoLP(:,:,c) = reshape(dolp, H, W);
        AoLP(:,:,c) = reshape(mod(aolp, pi), H, W);
    end
        plot_pol_parameters_comparison(a0, a1, a2, DoLP, AoLP, ' ', 'Fourier 4th');
end
