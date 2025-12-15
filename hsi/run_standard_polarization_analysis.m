function run_standard_polarization_analysis(dataTable, basePath, outputFolderOG) % Handles the Standard Stokes-based 4-angle dataset (WG @ 0°, 45°, 90°, 135°)

    % Select sample type (and day, when relevant)
    [filteredData, selectedType, selectedDay] = filter_data_by_sample_and_day(dataTable);

    disp('Computing reflectances for 4-angle dataset...');

    %% Compute mean reflectance per polarization angle
    % Sort by angle
    anglesDeg = cellfun(@str2double, filteredData.WGpolAngle);
    [sortedAngles, sortIdx] = sort(anglesDeg);
    filteredData = filteredData(sortIdx, :);

    mean_ref_struct = struct();
    spectral_ref_struct = struct();
    for i = 1:numel(sortedAngles)  %for each WG position/angle
        fpath = fullfile(basePath, filteredData.FolderName{i}, 'capture\');
        fname = filteredData.FolderName{i};

        % Extract the calibrated hyperspectral data (HS_calibrated) and region of interest (ROI)
        [Data, White, Dark, wavelengths] = read_data(fpath, fname);
        HS_calib = apply_calibration(Data, White, Dark);
        hsfiltered = apply_sg_filter(HS_calib);
        % Compute mean reflectance over the wavelengths
        mean_ref = mean(hsfiltered, 3); %x,y,lambda

        % Store results
        mean_ref_struct.(sprintf('angle_%d', sortedAngles(i))) = mean_ref;
        spectral_ref_struct.(sprintf('angle_%d', sortedAngles(i))) = single(hsfiltered); % reduce size
    end

    % Ensure same size for all angles
    mean_ref_struct = resize_reflectances(mean_ref_struct);
    cube = resize_hs_cubes(spectral_ref_struct); % also average across angles to reduce memory occupied

    outFolder = fullfile(outputFolderOG, [selectedType '_' selectedDay]);
    if ~exist(outFolder, 'dir'); mkdir(outFolder); end

    save_analysis_data(struct('Spectral_reflectances', cube,'Spatial_reflectances', mean_ref_struct, 'wavelengths', wavelengths), 'Standard_Reflectances', [selectedType '_' selectedDay], outFolder);
    
    %% Compute Stokes Parameters
    S0 = mean_ref_struct.angle_0 + mean_ref_struct.angle_90; % S0 = L0 + L90 (total intensity)
    S1 = mean_ref_struct.angle_0 - mean_ref_struct.angle_90; % S1 = L0 - L90 (horizontal – vertical)
    S2 = mean_ref_struct.angle_45 - mean_ref_struct.angle_135; % S2 = L45 - L135

    %% Standard computation of DoLP and AoLP as spatial functions 
    %DoLP = compute_dolp(mean_ref_struct);
    DoLP = sqrt(S1.^2 + S2.^2) ./ S0; % Compute Degree of Linear Polarization (DoLP)

    %AoLP = compute_aolp(mean_ref_struct);
    AoLP = 0.5 * atan2(S2, S1); % Compute Angle of Linear Polarization (AoLP) % atan2 for correct quadrant handling
    AoLP = rad2deg(AoLP); % Convert radians to degrees

    plot_polarization_spatially(DoLP, 'jet', ['DoLP (Standard) - ' selectedType], outFolder);
    plot_polarization_spatially(AoLP, 'hsv', ['AoLP (Standard) - ' selectedType], outFolder);

    %% Save results
    %save_analysis_data(struct('DoLP', DoLP, 'AoLP', AoLP), 'Standard_Stokes_Results', selectedType, outFolder);

    save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'DoLP', DoLP, 'AoLP', AoLP), 'Standard_Stokes_Results', selectedType, outFolder);

    disp(['Standard polarization analysis complete for ' selectedType]);
end
