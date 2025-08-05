function [S0, S1, S2, DoLP, AoP_deg] = compute_spie(I_theta)
    % Constants from calibration (from SPIE paper)
    C1 = 198.4;  % Radiometric calibration factor - value used for SPIE system
    % The camera has no intrinsic polarization sensitivity
    m12 = 0; % Horizontal polarimetric sensitivity
    m13 = 0; % Vertical polarimetric sensitivity
    
    % Flatten image cube for vectorized FFT (pixels x angles)
    [rows, cols, numAngles] = size(I_theta);
    Npix = rows * cols;
    I_pixels = reshape(I_theta, Npix, numAngles);
    
    % FFT across angles
    F = fft(I_pixels, [], 2) / numAngles;  % Normalized FFT
    
    DC = real(F(:,1));           % f=0 component
    F2 = F(:,3);                 % f=2: polarization signal
    F4 = F(:,5);                 % f=4: system response
    
    % Intrinsic offset & internal source known or estimated elsewhere
    % For simplicity, could be set to zero if unknown
    Si_int = 0; 
    offset = 0;
    
    % Build system of equations (Eq. 19 in SPIE paper)
    Y = [DC - (C1/2)*Si_int - offset,  real(F2) + (C1/2)*m12*Si_int, imag(F2) + (C1/2)*m13*Si_int, real(F4), imag(F4)];
    
    A = (C1/2) * [1,  m12,  m13; m12, 1, 0; m13, 0, 1; 0, m12/2, -m13/2; 0, m13/2,  m12/2];
    
    % Solve for S0, S1, S2 (vectorized across pixels)
    S = pinv(A) * Y';  % [3 x Npix]
    
    % Reshape to 2D images
    S0 = reshape(S(1,:), rows, cols);
    S1 = reshape(S(2,:), rows, cols);
    S2 = reshape(S(3,:), rows, cols);
    
    % Compute Polarization Parameters
    DoLP = sqrt(S1.^2 + S2.^2) ./ S0;
    AoP_rad = 0.5 * atan2(S2, S1);
    AoP_deg = rad2deg(AoP_rad);
   
end
