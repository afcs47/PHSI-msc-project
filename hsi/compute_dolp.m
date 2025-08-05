function DoLP = compute_dolp(mean_reflectances)
    % Compute Stokes parameters
    S0 = mean_reflectances.('angle_0') + mean_reflectances.('angle_90'); % S0 = L0 + L90
    S1 = mean_reflectances.('angle_0') - mean_reflectances.('angle_90'); % S1 = L0 - L90
    S2 = mean_reflectances.('angle_45') - mean_reflectances.('angle_135'); % S2 = L45 - L135
    
    % Compute Degree of Linear Polarization (DoLP)
    DoLP = sqrt(S1.^2 + S2.^2) ./ S0; 
end