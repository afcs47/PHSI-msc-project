function [DoLP, AoLP] = calculate_polarization_fourier_2nd(polar_images, datasetname, outputFolder) % outputs DoLP, AoLP, computed from 2nd order Fourier model from the inputed polar_images demosaiced polarization stack
    % Define angles used in the camera
    angles_deg = [0, 45, 90, 135];

    % Create the matrix for  2nd-order fit
    A2 = [ones(4,1), cosd(2 * angles_deg)', sind(2 * angles_deg)']; %  2nd-order model: [a0, a2, b2]

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

        % Solve linear least squares for 2nd-order coefficients
        coeffs = A2 \ I_reshaped; % [3, H*W]
        a0 = coeffs(1, :);
        a2 = coeffs(2, :);
        b2 = coeffs(3, :);

         % Compute DoLP and AoLP from Fourier coefficients
        dolp = sqrt(a2.^2 + b2.^2) ./ (a0 + eps); % avoid divide by zero
        aolp = 0.5 * atan2(b2, a2); % angle in radians
        
       % Reshape back to [H, W]
        DoLP(:,:,c) = reshape(dolp, H, W);
        AoLP(:,:,c) = reshape(aolp, H, W); % still in radians

    end

    plot_pol_parameters_comparison(reshape(a0, H, W), reshape(a2, H, W), reshape(b2, H, W), DoLP, AoLP, datasetname, 'Fourier 2nd', outputFolder);
end
