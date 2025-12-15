function [hsi_rgb_u8, pol_rgb] = prepare_calibration_images(basePath, calibHsData, calibPolData, selectedType, selectedDay, datasetname)

    % Select calibration HSI and pol data 
    selectedCalibData = calibHsData(strcmp(calibHsData.SampleName, selectedType) & strcmp(calibHsData.Day, selectedDay), :);
    fname = selectedCalibData.FolderName{1};
    fpath = fullfile(basePath, fname, 'capture\');
    
    names = calibPolData.FileName; %2×1 cell array, since pol was measured before and after the hsi measurement for reference
    badPattern = selectedType + "2_"; % build the unwanted pattern (eg: solutions2_)
    filteredName = names(~contains(names, badPattern)& contains(names, selectedDay)); % keep only the entries that do mot contain that bad pattern

    calib_pol_raw = fullfile(basePath, "pol data/", filteredName);


    %% Load and calibrate HSI image 
    [Data, White, Dark, ~] = read_data(fpath, fname);
    HS_calibrated = apply_calibration(Data, White, Dark);

    [R, G, B] = fake_rgb(false, HS_calibrated, datasetname); % Get fake RGB band indices
    hsfiltered = apply_sg_filter(HS_calibrated); % Spectral smoothing
    hsi_rgb = hsfiltered(:,:,[R, G, B]); % Compose the pseudo-RGB image
    hsi_rgb_u8 = im2uint8(hsi_rgb); % Convert to uint8 for display

    %% Load and process the polarization image
    row = 2048; col = 2448;
    fin = fopen(calib_pol_raw, 'r');
    I = fread(fin, row * col, 'uint8=>uint8');
    Z = reshape(I, col, row)';
    fclose(fin);
    %J = demosaic(Z, "rggb"); % Demosaic raw Bayer image

    fun = @(block) mean(mean(block.data));
    proc = blockproc(Z, [2 2], fun); % Average 2x2 blocks to reduce artifacts
    pol_rgb = demosaic(uint8(proc), 'rggb');
    pol_rgb = imrotate(pol_rgb, 180); % Correct the 180 deg rotation
    %pol_rgb = rot90(pol_rgb, 2);  % Equivalent to 180 deg rotation

end

