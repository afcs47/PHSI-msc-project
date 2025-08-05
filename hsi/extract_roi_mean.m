function [mean_ref_fc, std_ref] = extract_roi_mean(I_roi)
    sz = size(I_roi); 
    I = reshape(I_roi, [sz(1) * sz(2), sz(3)]); % Reshapes data for processing
    mean_ref_fc = mean(I); % Computes the mean reflectance
    std_ref = std(I); % Computes the standard deviation
end