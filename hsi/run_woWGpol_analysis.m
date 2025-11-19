function run_woWGpol_analysis(dataTable, basePath, outputFolderOG) % Handles datasets without waveplate (woWGpolData) ie computes reflectance-only analysis

    %Select sample type (and day, when relevant)
    [filteredData, selectedType, selectedDay] = filterDataBySampleAndDay(dataTable);

    disp('Computing reflectance data for no-WG polarizer dataset...');

    %% Compute mean reflectance
    mean_ref_struct = struct();
    spectral_ref_struct = struct();
    for i = 1:height(filteredData)
        fpath = fullfile(basePath, filteredData.FolderName{i}, 'capture\');
        fname = filteredData.FolderName{i};

        % Extract the calibrated hyperspectral data
        [Data, White, Dark, wavelengths] = read_data(fpath, fname);
        HS_calib = apply_calibration(Data, White, Dark);
        hsfiltered = apply_sg_filter(HS_calib);
        % Compute mean reflectance over the wavelengths
        mean_ref = mean(hsfiltered, 3);

        % Store results
        mean_ref_struct.(sprintf('sample_%02d', i)) = mean_ref;
        spectral_ref_struct.(sprintf('sample_%02d', i)) = hsfiltered;
    end

    %% Plot reflectance
    plot_reflectances_spatially(mean_ref, 'gray', ['Normalized Reflectance for ' selectedType selectedDay 'junwoWGpol (grayscale)' ], fullfile(outputFolderOG, [selectedType selectedDay 'junwoWGpol']));
    plot_reflectances_spatially(mean_ref, 'jet', ['Normalized Reflectance for ' selectedType selectedDay 'junwoWGpol'], fullfile(outputFolderOG, [selectedType selectedDay 'junwoWGpol']));

    %% Save results
    outFolder = fullfile(outputFolderOG, [selectedType '_' selectedDay '_woWG']);
    if ~exist(outFolder, 'dir'); mkdir(outFolder); end

    save_analysis_data(struct('Spectral_reflectances', spectral_ref_struct, 'Spatial_reflectances', mean_ref_struct, 'wavelengths', wavelengths), 'woWG_Reflectance_Results', selectedType, outFolder); %reflectances over wavelengths and over spatial coordinates

    disp(['woWGpolData analysis complete and saved in: ' outFolder]);
end
