function run_standard_polarization_analysis(dataTable, basePath, outputFolderOG) % Handles the Standard Stokes-based 4-angle dataset (WG @ 0°, 45°, 90°, 135°)

    % Select sample type (and day, when relevant)
    [filteredData, selectedType, selectedDay] = filterDataBySampleAndDay(dataTable);

    disp('Computing reflectances for 4-angle dataset...');

    %% Compute mean reflectance per polarization angle
    anglesDeg = cellfun(@str2double, filteredData.WGpolAngle);
    [sortedAngles, sortIdx] = sort(anglesDeg);
    filteredData = filteredData(sortIdx, :);

    mean_ref_struct = struct();
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
    end

    % Ensure same size for all angles
    mean_ref_struct = resize_reflectances(mean_ref_struct);

    outFolder = fullfile(outputFolderOG, [selectedType '_' selectedDay]);
    if ~exist(outFolder, 'dir'); mkdir(outFolder); end

    %% Standard computation of DoLP and AoLP as spatial functions 
    DoLP = compute_dolp(mean_ref_struct);
    AoLP = compute_aolp(mean_ref_struct);

    plot_polarization_spatially(DoLP, 'jet', ['DoLP (Standard) - ' selectedType], outFolder);
    plot_polarization_spatially(AoLP, 'hsv', ['AoLP (Standard) - ' selectedType], outFolder);

    %% Save results
    save_analysis_data(struct('DoLP', DoLP, 'AoLP', AoLP), 'Standard_Stokes_Results', selectedType, outFolder);

    disp(['Standard polarization analysis complete for ' selectedType]);
end
