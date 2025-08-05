function [S0, S1, S2, DoLP_fourier, AoLP_fourier] = compute_fourier_pol(I_theta, sortedAngles, numAngles)
    % Convert angle degrees to radians and prepare cosine/sine vectors
    theta_rad = deg2rad(sortedAngles(:)); % Column vector
    cos2theta = cos(2 * theta_rad);       % [numAngles x 1]
    sin2theta = sin(2 * theta_rad);
    
    % Reshape intensity data into [pixels x angles] for vectorized processing
    I_pixels = reshape(I_theta, [], numAngles); % [numPixels x numAngles]
    
    % Fourier fit: compute a0, a1, a2 via least square method
    A = [ones(numAngles,1), cos2theta, sin2theta]; % Design matrix [numAngles x 3]
    coeffs = A \ I_pixels.'; % Solve: [3 x numPixels] (transposed for better conditioning)
    
    a0 = coeffs(1, :).'; % [numPixels x 1]
    a1 = coeffs(2, :).';
    a2 = coeffs(3, :).';
    
    % Compute Stokes parameters
    S0 = 2 * a0; % Total intensity
    S1 = 2 * a1;
    S2 = 2 * a2;
    
    % Compute DoLP and AoLP
    DoLP_fourier = sqrt(S1.^2 + S2.^2) ./ S0; % [numPixels x 1]
    AoLP_fourier = 0.5 * atan2(S2, S1); % Radians

end