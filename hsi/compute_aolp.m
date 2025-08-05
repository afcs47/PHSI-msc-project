
function AoLP = compute_aolp(mean_reflectances)
    % Compute Stokes parameters
    S1 = mean_reflectances.('angle_0') - mean_reflectances.('angle_90'); % S1 = L0 - L90
    S2 = mean_reflectances.('angle_45') - mean_reflectances.( 'angle_135'); % S2 = L45 - L135

    % Compute Angle of Linear Polarization (AoLP)
    AoLP = 0.5 * atan2(S2, S1); % atan2 for correct quadrant handling

    % Convert radians to degrees
    AoLP = rad2deg(AoLP);
end