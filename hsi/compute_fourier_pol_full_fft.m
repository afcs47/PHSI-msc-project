function [S0, S1, S2, S3, DoLP_fft, AoLP_fft, Ellipticity] = compute_fourier_pol_full_fft(I_theta, numAngles)
    % Ensure data is [pixels x numAngles]
    if size(I_theta,2) ~= numAngles
        I_pixels = reshape(I_theta, [], numAngles);
    else
        I_pixels = I_theta;
    end
    
    % Apply FFT along angle dimension
    F = fft(I_pixels, [], 2) / numAngles; % normalize
    
    % Extract Fourier coefficients:
    % F(:,1) - DC term (a0)
    % F(:,2) - harmonic at 1 cycle (cos2θ/sin2θ contributions)
    % F(:,3) - harmonic at 2 cycles (cos4θ/sin4θ contributions)
    
    a0 = 2*real(F(:,1)); % DC component
    c2 = 2*real(F(:,2)); % cos(2θ)
    s2 = -2*imag(F(:,2)); % sin(2θ)
    s4 = -2*imag(F(:,3)); % sin(4θ) (linked to S3)
    
    % Map Fourier coefficients to Stokes parameters
    S0 = a0;
    S1 = c2;
    S2 = s2;
    S3 = s4;
    
    % Derived polarization quantities
    DoLP_fft = sqrt(S1.^2 + S2.^2) ./ S0;
    AoLP_fft = 0.5 * atan2(S2, S1);
    Ellipticity = asin(S3 ./ S0);
end
