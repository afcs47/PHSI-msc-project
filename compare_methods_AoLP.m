function compare_methods_AoLP(mean_AoLP_standard, method1, mean_AoLP_1, method2, mean_AoLP_2, method3, min_angle, max_angle)
    
%     % Convert range from degrees to radians
     %min_rad = deg2rad(min_angle);
     %max_rad = deg2rad(max_angle);
%     range_rad = max_rad - min_rad;
% 
%     % Wrap all AoLP maps into [min_rad, max_rad)
%     wrap_range = @(x) mod(x - min_rad, range_rad) + min_rad;
%     mean_AoLP_standard = wrap_range(mean_AoLP_standard);
%     mean_AoLP_1 = wrap_range(mean_AoLP_1);
%     mean_AoLP_2 = wrap_range(mean_AoLP_2);
%    

    % Convert each to HSV-based RGB image to use the same color scheme
    rgb_standard = map_aolp_to_rgb(mean_AoLP_standard, min_angle, max_angle);
    % rgb_1 = map_aolp_to_rgb(mean_AoLP_1, min_angle, max_angle);
    % rgb_2 = map_aolp_to_rgb(mean_AoLP_2, min_angle, max_angle);
    
    % figure('Name', 'AoLP Comparison (HSV Hue Mapping)');
    % subplot(1,3,1);
    % imshow(rgb_standard);
    % title(sprintf('AoLP %s (deg)', method1));
    % 
    % subplot(1,3,2);
    % imshow(rgb_1);
    % title(sprintf('AoLP %s (deg)', method2));
    % 
    % subplot(1,3,3);
    % imshow(rgb_2);
    % title(sprintf('AoLP %s (deg)', method3));
    
    % Compute angular differences (in radians) - Shortest Angular Path (AoLP is a circular variable (angle), not linear / angle difference must be circular ->e.g.: 179deg vs 1deg is a 2deg diff, not 178deg)
    diff_1 = atan2(sin(mean_AoLP_1 - mean_AoLP_standard), cos(mean_AoLP_1 - mean_AoLP_standard));
    diff_2 = atan2(sin(mean_AoLP_2  - mean_AoLP_standard), cos(mean_AoLP_2 - mean_AoLP_standard));
    diff_fourier = atan2(sin(mean_AoLP_1 - mean_AoLP_2 ), cos(mean_AoLP_1 - mean_AoLP_2));
    
    % Convert to degrees
    diff_1_deg = rad2deg(diff_1); % Range: [-90deg, +90deg]
    diff_2_deg = rad2deg(diff_2);
    diff_fourier_deg = rad2deg(diff_fourier);
    
    figure('Name', 'AoLP Comparison');
    
    subplot(2,2,1);
    imshow(rgb_standard);
    colormap(hsv);
    set_aolp_colorbar(gca, min_angle, max_angle);
    title(sprintf('%s AoLP (deg)', method1)); 
    
    subplot(2,2,2);
    %imagesc(diff_1_deg);
    imshow(diff_1_deg, []);
    colormap(gca, 'hsv'); colorbar;
    set_aolp_colorbar(gca, min_angle, max_angle);
    title(sprintf('%s - %s (deg)', method2, method1));
    
    subplot(2,2,3);
    %imagesc(diff_2_deg);
    imshow(diff_2_deg, []);
    colormap(gca, 'hsv'); colorbar;
    set_aolp_colorbar(gca, min_angle, max_angle);
    title(sprintf('%s - %s (deg)', method3, method1));
    
    
    subplot(2,2,4);
    %imagesc(diff_fourier_deg);
    imshow(diff_fourier_deg, []);
    colormap(gca, 'hsv'); colorbar;
    set_aolp_colorbar(gca, min_angle, max_angle);
    title(sprintf('%s - %s (deg)', method2, method3));

end