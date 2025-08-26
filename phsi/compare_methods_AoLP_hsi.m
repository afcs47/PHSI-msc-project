function compare_methods_AoLP_hsi(mean_AoLP_standard, method1, mean_AoLP_1, method2, mean_AoLP_2, method3)
    figure('Name', 'AoLP Comparison');
    
    subplot(2,2,1);
    imshow(mean_AoLP_standard, []); 
    colormap(hsv); colorbar;
    title(sprintf('%s AoLP (deg)', method1)); 

    % Wrap AoLP to [0, pi] before comparison
%     mean_AoLP_standard = mod(mean_AoLP_standard, pi);
%     mean_AoLP_1= mod(mean_AoLP_fourier1, pi);
%     mean_AoLP_2 = mod(mean_AoLP_fourier2, pi);
    
    % Compute angular differences (in radians) - Shortest Angular Path (AoLP is a circular variable (angle), not linear / angle difference must be circular ->e.g.: 179deg vs 1deg is a 2deg diff, not 178deg)
    diff_1 = atan2(sin(mean_AoLP_1 - mean_AoLP_standard), cos(mean_AoLP_1 - mean_AoLP_standard));
    diff_2 = atan2(sin(mean_AoLP_2  - mean_AoLP_standard), cos(mean_AoLP_2 - mean_AoLP_standard));
    diff_fourier = atan2(sin(mean_AoLP_1 - mean_AoLP_2 ), cos(mean_AoLP_1 - mean_AoLP_2));
    
    % Convert to degrees
    diff_1_deg = rad2deg(diff_1); % Range: [-90deg, +90deg]
    diff_2_deg = rad2deg(diff_2);
    diff_fourier_deg = rad2deg(diff_fourier);
    
    
    subplot(2,2,2);
    %imagesc(diff_1_deg);
    imshow(diff_1_deg, []);
    colormap(gca, 'hsv'); colorbar;
    caxis([-90, 90]);
    title(sprintf('%s - %s (deg)', method2, method1));
    
    subplot(2,2,3);
    %imagesc(diff_2_deg);
    imshow(diff_2_deg, []);
    colormap(gca, 'hsv'); colorbar;
    caxis([-90, 90]);
    title(sprintf('%s - %s (deg)', method3, method1));
    
    
    subplot(2,2,4);
    %imagesc(diff_fourier_deg);
    imshow(diff_fourier_deg, []);
    colormap(gca, 'hsv'); colorbar;
    caxis([-90, 90]);
    title(sprintf('%s - %s (deg)', method2, method3));

end