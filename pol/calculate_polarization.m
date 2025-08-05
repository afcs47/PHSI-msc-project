%{
calculate_polarization - is used as a helper function to calculate DoLP and
AoLP, to reduce artifacts in the calculated image we demosaic each set of
polarized pixels separately. Then we calculate the Stokes vector from given
intensities and we use those to calculate DoLP and AoLP for each pixel

@Z - raw image to use for processing
%}

function [DoLP, AoLP] = calculate_polarization(pol_90, pol_45, pol_135, pol_0) % Standard Stokes parameters calculation - switch input variables to individual polarization images to avoid repeated computation; DoLP computed in an alternative way 
    % demosaic each polarization separately
    color_90 = demosaic(uint8(pol_90), 'rggb');
    color_45 = demosaic(uint8(pol_45), 'rggb');
    color_135 = demosaic(uint8(pol_135), 'rggb');
    color_0 = demosaic(uint8(pol_0), 'rggb');
    
    [height, width, channels] = size(color_0);
    
    % init DoLP and AoLP
    DoLP = zeros(height, width, channels);
    AoLP = zeros(height, width, channels);
    
    % process each color channel
    for c = 1:channels
        % intensities for each polarization angle (need to be type double to
        % work in sqrt function)
        I_0 = double(color_0(:,:,c));
        I_90 = double(color_90(:,:,c));
        I_45 = double(color_45(:,:,c));
        I_135 = double(color_135(:,:,c));
        
        %Stokes vector
        S0 = I_0 + I_90;  % total image intensity
        S1 = I_0 - I_90;  % diff between horizontal and vertical polarization
        S2 = I_45 - I_135; % diff between 45° and 135° polarization
        
        %DoLP(:,:,c) = sqrt(S1.^2 + S2.^2) ./ (S0 + eps); %eps to not divide by zero

        %Alternative DoLP computation: Uses a threshold mask on S0 to
        %compute DoLP only where the total intensity is reliable, avoiding
        %extreme values in dark/noisy pixels -> the extreme values go from
        %[0.00 36028797018963968.00] to [0.00 2.00] without affecting the
        %resulting DoLP representation
        valid = S0 > 1;
        DoLP_channel = zeros(size(S0));
        DoLP_channel(valid) = sqrt(S1(valid).^2 + S2(valid).^2) ./ S0(valid);
        DoLP(:,:,c) = DoLP_channel;
        
        AoLP(:,:,c) = 0.5 * atan2(S2, S1);
    end

    plot_pol_parameters_comparison(S0, S1, S2, DoLP, AoLP, ' ', 'Standard');
end