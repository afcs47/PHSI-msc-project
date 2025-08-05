function HS_calibrated = apply_calibration(Data, White, Dark)
    HS_calibrated = zeros(size(Data)); % Initializes a zero matrix for the calibrated data
    white_ref = mean(White, 1); % Averages the white reference across the first dimension
    dark_ref = mean(Dark, 1); % Averages the dark reference across the first dimension
    
    % Calibrates the hyperspectral data by applying the standard formula
    for i = 1:size(Data, 1)
       HS_calibrated(i,:,:) = (Data(i,:,:) - dark_ref(1,:,:))./(white_ref(1,:,:) - dark_ref(1,:,:));
    end
end