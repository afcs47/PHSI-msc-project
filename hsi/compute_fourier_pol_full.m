function [S0, S1, S2, S3, DoLP_fourier, AoLP_fourier, Ellipticity] = compute_fourier_pol_full(I_theta, sortedAngles, numAngles)
    % Convert angle degrees to radians
    theta_rad = deg2rad(sortedAngles(:)); % [numAngles x 1]
    
    % Build Fourier basis up to 4θ harmonics
    cos2theta = cos(2 * theta_rad);
    sin2theta = sin(2 * theta_rad);
    cos4theta = cos(4 * theta_rad);
    sin4theta = sin(4 * theta_rad);
    
    % Reshape intensity data into [pixels x angles]
    I_pixels = reshape(I_theta, [], numAngles); % [numPixels x numAngles]
    
    % Design matrix with harmonics (DC, cos2θ, sin2θ, cos4θ, sin4θ)
    A = [ones(numAngles,1), cos2theta, sin2theta, cos4theta, sin4theta];
    
    % Solve least-squares system
    coeffs = A \ I_pixels.'; % [5 x numPixels]
    
    a0 = coeffs(1, :).'; % DC term
    a1 = coeffs(2, :).'; % cos(2θ) - S1
    a2 = coeffs(3, :).'; % sin(2θ) - S2
    a3 = coeffs(4, :).'; % cos(4θ) - contributes to S1/S2 calibration
    a4 = coeffs(5, :).'; % sin(4θ) - linked to S3
    
    % Map Fourier coefficients to Stokes parameters
    % Note: scaling depends on convention; here we follow Berry & Schaefer
    S0 = 2 * a0;
    S1 = 2 * a1;
    S2 = 2 * a2;
    S3 = 2 * a4; % circular polarization from sin(4θ) harmonic
    
    % Derived polarization parameters
    DoLP_fourier = sqrt(S1.^2 + S2.^2) ./ S0;
    AoLP_fourier = 0.5 * atan2(S2, S1); % radians
    Ellipticity  = asin(S3 ./ S0); % ellipticity angle (optional)
end
