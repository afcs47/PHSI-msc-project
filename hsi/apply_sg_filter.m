function hsfiltered = apply_sg_filter(HS_calibrated)
    hsfiltered = zeros(size(HS_calibrated));  % Initializes the filtered data matrix
    order = 2; % Set order for SG filter
    window = 15; % Set size of window of spectral channels for SG filter 
    
    % Apply the filter to each pixel's spectral data
    for x = 1:length(HS_calibrated(:,1,1))
        for y = 1:length(HS_calibrated(1,:,1))
            hsfiltered(x, y, :) = sgolayfilt(HS_calibrated(x, y, :), order, window);
        end
    end
end