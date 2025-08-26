function tform = calibrateHsiPolData(fpath, fname_selected, pol_raw_path, datasetname, resizeSizeHsi, resizeSizePol, check_calibration) %Aligns and fuses HSI and Polarization images via checkerboard homography
    % Add auxiliary HSI functions to path
    addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi'));
    
    %% Load and calibrate HSI image
    [Data, White, Dark, wavelengths] = read_data(fpath, fname_selected);
    HS_calibrated = apply_calibration(Data, White, Dark);
    [R, G, B] = fake_rgb(false, HS_calibrated, datasetname); % Get fake RGB band indices
    hsfiltered = apply_sg_filter(HS_calibrated); % Spectral smoothing
    hsi_rgb = hsfiltered(:,:,[R, G, B]); % Compose the pseudo-RGB image
    hsi_rgb = hsi_rgb(1:resizeSizeHsi(1), 1:resizeSizeHsi(2), :); % Resize to the size of DoLP/AoLP to which apply tform 
    hsi_rgb_u8 = im2uint8(hsi_rgb); % Convert to uint8 for display
    
    %% Load and process the polarization image
    row = 2048; col = 2448;
    fin = fopen(pol_raw_path, 'r');
    I = fread(fin, row * col, 'uint8=>uint8');
    Z = reshape(I, col, row)';
    fclose(fin);
    %J = demosaic(Z, "rggb"); % Demosaic raw Bayer image
    fun = @(block) mean(mean(block.data));
    proc = blockproc(Z, [2 2], fun); % Average 2x2 blocks to reduce artifacts
    pol_rgb = demosaic(uint8(proc), 'rggb');
    pol_rgb = imrotate(pol_rgb, 180); % Correct the 180 deg rotation
    %pol_rgb = rot90(pol_rgb, 2);  % Equivalent to 180 deg rotation
    pol_rgb = pol_rgb(1:resizeSizePol(1), 1:resizeSizePol(2), :); % Resize to the size of DoLP/AoLP to which apply tform 
    
    %% Attempt checkerboard detection
    fprintf("Attempting automatic checkerboard detection...\n");
    % Try multiple enhancements for checkerboard visibility
    [imagePointsHSI, boardSizeHSI] = tryDetectCheckerboard(hsi_rgb_u8, 'HSI'); % Detect in HSI RGB
    [imagePointsPOL, boardSizePOL] = tryDetectCheckerboard(pol_rgb, 'Polarization'); % Detect in Polarization RGB
    
    if isempty(imagePointsHSI) || isempty(imagePointsPOL)
        error("Checkerboard could not be reliably detected in one or both images.");
    end
    
    % Display detected corners
    % figure, imshow(hsi_rgb_u8), hold on;
    % plot(imagePointsHSI(:,1), imagePointsHSI(:,2), 'ro'), title('HSI Checkerboard Corners');
    % figure, imshow(pol_rgb), hold on;
    % plot(imagePointsPOL(:,1), imagePointsPOL(:,2), 'go'), title('Pol Checkerboard Corners');
    
    %% Select 4 matching points by clicking on the detected green dots
    nPointsToSelect = 4;
    
    % HSI Image
    figure; imshow(hsi_rgb_u8); axis image; % Ensure square pixels
    title('HSI: Click 4 checkerboard points (close to the green marks)');
    hold on; plot(imagePointsHSI(:,1), imagePointsHSI(:,2), 'g+');
    text(imagePointsHSI(:,1)+3, imagePointsHSI(:,2), arrayfun(@num2str, 1:size(imagePointsHSI,1), 'UniformOutput', false), 'Color', 'green', 'FontSize', 8);
    hold off;
    
    ptsHSI = zeros(nPointsToSelect, 2);
    for i = 1:nPointsToSelect
        [x, y] = ginput(1);
        % Find nearest detected point to get a better corner selection accuracy and precision
        dists = vecnorm(imagePointsHSI - [x y], 2, 2);
        [~, idx] = min(dists);
        ptsHSI(i, :) = imagePointsHSI(idx, :);
        % Show selection
        hold on; plot(imagePointsHSI(idx,1), imagePointsHSI(idx,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    % Polarization Image
    figure; imshow(pol_rgb); axis image; % Ensure square pixels
    title('POL: Click the same 4 checkerboard points (in the same order)');
    hold on; plot(imagePointsPOL(:,1), imagePointsPOL(:,2), 'g+');
    text(imagePointsPOL(:,1)+3, imagePointsPOL(:,2), arrayfun(@num2str, 1:size(imagePointsPOL,1), 'UniformOutput', false), 'Color', 'green', 'FontSize', 8);
    hold off;
    
    ptsPOL = zeros(nPointsToSelect, 2);
    for i = 1:nPointsToSelect
        [x, y] = ginput(1);
        % Find nearest detected point
        dists = vecnorm(imagePointsPOL - [x y], 2, 2);
        [~, idx] = min(dists);
        ptsPOL(i, :) = imagePointsPOL(idx, :);
        % Show selection
        hold on; plot(imagePointsPOL(idx,1), imagePointsPOL(idx,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    %% Estimate homography and warp polarization image

    % Since there are 2 possible Matlab functions available to use, and to 
    % make sure the misalignment issue doesn't come from the homography 
    % matrix used, both methods will be applied, with the best one being chosen at the end

    % estgeotform2d - Robust, with outlier rejection, good for autodetected points
    % fitgeotform2d - Direct least-squares fit, good for 4 manually selected points
    
    [tform_best, tform_fit,tform_est, report] = compareCalibrationMethodsROI(ptsHSI, ptsPOL, imagePointsHSI, imagePointsPOL,hsi_rgb_u8, pol_rgb, check_calibration);

    fprintf('\nRMSE fitgeotform2d = %.3f px, estgeotform2d = %.3f px\n\n', report.err_fit, report.err_est);
    
    [bestErr, bestIdx] = min([report.err_fit, report.err_est]);
    names = {'tform_fit','tform_est'};
    % Menu selection
    choice = menu(sprintf('Select transformation (recommended: %s)', names{bestIdx}), 'Use tform_fit','Use tform_est');
    % Assign chosen transform
    tform = eval(names{choice});  
    fprintf('Using %s.\n', names{choice});

    % Use the best transform for everything else:
    alignedPOL = imwarp(pol_rgb, tform, 'OutputView', imref2d(size(hsi_rgb_u8)));

    %% Show and save fused result (tranform matrix) for future use
    figure; imshowpair(hsi_rgb_u8, alignedPOL, 'falsecolor'); %Colored regions show where the intensities differ
    title('Fused Image: HSI + aligned POL');
    axis image; % Ensure square pixels
    
    save([datasetname '_homography_tform.mat'], 'tform');
    disp("✅ Calibration completed. Parameters saved to 'homography_tform.mat'."); 

    if check_calibration
        disp('For the detected checkerboard points:');
        imagePointsHSI 
        boardSizeHSI
        imagePointsPOL
        boardSizePOL

        disp('The following corners were selected:')
        disp('HSI points:'), disp(ptsHSI);
        disp('POL points:'), disp(ptsPOL);

        figure; imshowpair(hsi_rgb_u8,alignedPOL,"diff");
        title('Difference between HSI and aligned POL');
        axis image; % Ensure square pixels

        %% Compute residuals from 2D point mapping
        % Transform POL points to HSI image using estimated homography
        ptsPOL_warped = transformPointsForward(tform, ptsPOL);
        
        % Compute 2D residuals
        residuals2D = vecnorm(ptsHSI - ptsPOL_warped, 2, 2) % Euclidean distance per point
        
        fprintf("Mean 2D reprojection error: %.4f px\n", mean(residuals2D));
        fprintf("STD of 2D reprojection error: %.4f px\n", std(residuals2D));

        maxDist = max(vecnorm(ptsHSI, 2, 2));
        residualsPercentOfImage = residuals2D / maxDist
        residuals2DToOrigin = residuals2D./vecnorm(ptsHSI,2,2)

        %% Evaluate Fusion Precision based on matched checkerboard corner points

        % Apply homographic transformation
        ptsPOL_proj = transformPointsForward(tform, ptsPOL);
        
        % Calculate Euclidean distance error
        distances = sqrt(sum((ptsHSI - ptsPOL_proj).^2, 2));
        fusion_precision = mean(distances);  % in pixels
        
        fprintf('Fusion Precision: %.4f pixels (average Euclidean distance)\n', fusion_precision);
        
        %% Transform all POL checkerboard points into the full HSI image space
        % Get number of matching points
        nPoints = min(size(imagePointsHSI, 1), size(imagePointsPOL, 1));
        
        % Truncate to the same number of points
        ptsHSI_full = imagePointsHSI(1:nPoints, :);
        ptsPOL_full = imagePointsPOL(1:nPoints, :);
        
        % Warp the POL points to HSI space
        ptsPOL_warped_full = transformPointsForward(tform, ptsPOL_full);
        
        % Compute residuals
        residualsFull = vecnorm(ptsHSI_full - ptsPOL_warped_full, 2, 2);
        
        maxRadialDist = max(vecnorm(imagePointsHSI, 2, 2));
        residualsPercent = residualsFull / maxRadialDist * 100;
        
        % Report basic stats
        fprintf("Full checkerboard reprojection error:");
        fprintf("\n  Mean = %.4f px (%.4f %%)", mean(residualsFull), mean(residualsPercent));
        fprintf("\n  STD  = %.4f px", std(residualsFull));
        fprintf("\n  Max  = %.4f px\n", max(residualsFull));
        
    end
end    

