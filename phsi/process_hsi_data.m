function [wavelengths, pos_rect, mean_ref, std_ref] = process_hsi_data(fpath, fname_selected, datasetname, showRGB, k, pos_rect)
    % Read HS image
    [Data, White, Dark, wavelengths] = read_data(fpath, fname_selected);
    
    % Calibration
    HS_calibrated = apply_calibration(Data, White, Dark);
    
    % Show fake RGB image
    [R, G, B] = fake_rgb(showRGB, HS_calibrated, datasetname);
    
    % Savitzky-Golay Filtering
    hsfiltered = apply_sg_filter(HS_calibrated);

    % Selection of region of interest (ROI)
    if k == 1 % First dataset determines the ROI
        validROI = false; % Flag to track if ROI is valid
        while ~validROI
            [validROI, pos_rect] = select_roi(hsfiltered,R,G,B);
        end
    end
    I_roi = hsfiltered(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)),:); % area of image

    % Extract the ROI data from the hs image and calculate mean reflectance from it
    [mean_ref, std_ref] = extract_roi_mean(I_roi);
end
