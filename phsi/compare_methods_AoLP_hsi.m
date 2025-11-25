function compare_methods_AoLP_hsi(mean_AoLP_standard, method1, mean_AoLP_1, method2, mean_AoLP_2, method3) %Compute angular differences between the different aolps 
%NOTE: Shortest Angular Path is used (AoLP is a circular variable (angle), not linear => angle difference must be circular -> e.g.: 179deg vs 1deg is a 2deg diff, not 178deg)
    figure('Name', 'AoLP Comparison');
    
    subplot(2,2,1);
    imshow(mean_AoLP_standard, []); 
    colormap(hsv); colorbar;
    title(sprintf('%s AoLP (deg)', method1)); 

    % Wrap AoLP to [0, pi] before comparison
%     mean_AoLP_standard = mod(mean_AoLP_standard, pi);
%     mean_AoLP_1= mod(mean_AoLP_fourier1, pi);
%     mean_AoLP_2 = mod(mean_AoLP_fourier2, pi);
    
%     % Compute angular differences (in radians) 
%     diff_1 = atan2(sin(mean_AoLP_1 - mean_AoLP_standard), cos(mean_AoLP_1 - mean_AoLP_standard));
%     diff_2 = atan2(sin(mean_AoLP_2  - mean_AoLP_standard), cos(mean_AoLP_2 - mean_AoLP_standard));
%     diff_12 = atan2(sin(mean_AoLP_1 - mean_AoLP_2 ), cos(mean_AoLP_1 - mean_AoLP_2));
%     
%     % Convert to degrees
%     diff_1_deg = rad2deg(diff_1); % Range: [-90deg, +90deg]
%     diff_2_deg = rad2deg(diff_2);
%     diff_12_deg = rad2deg(diff_12);
%     

    % Convert AoLP maps (ALREADY IN DEGREES) to radians before the comparison
    A_std_rad = deg2rad(mean_AoLP_standard);
    A_1_rad   = deg2rad(mean_AoLP_1);
    A_2_rad   = deg2rad(mean_AoLP_2);

    % Compute circular angular differences (radians)
    diff_1  = atan2(sin(A_1_rad - A_std_rad), cos(A_1_rad - A_std_rad) );
    diff_2  = atan2(sin(A_2_rad - A_std_rad), cos(A_2_rad - A_std_rad) );
    diff_12 = atan2(sin(A_1_rad - A_2_rad), cos(A_1_rad - A_2_rad) );

    % Convert result back to degrees
    diff_1_deg  = rad2deg(diff_1);
    diff_2_deg  = rad2deg(diff_2);
    diff_12_deg = rad2deg(diff_12);
    
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
    %imagesc(diff_12_deg);
    imshow(diff_12_deg, []);
    colormap(gca, 'hsv'); colorbar;
    caxis([-90, 90]);
    title(sprintf('%s - %s (deg)', method2, method3));

end