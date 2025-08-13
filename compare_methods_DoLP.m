function compare_methods_DoLP(mean_DoLP_standard, method1, mean_DoLP_fourier1, method2, mean_DoLP_fourier2, method3)
    diff_DoLP_1 = abs(mean_DoLP_standard - mean_DoLP_fourier1);
    diff_DoLP_2 = abs(mean_DoLP_standard - mean_DoLP_fourier2);
    diff_DoLP_fourier = abs(mean_DoLP_fourier1 - mean_DoLP_fourier2);
    
    figure('Name', 'DoLP Comparison');
    subplot(2,2,1);
    imshow(mean_DoLP_standard, []);
    colormap(jet); colorbar;
    caxis([0 1]);  % or set(gca, 'CLim', [0 1])
    title(sprintf('%s DoLP', method1));   
    
    subplot(2,2,2);
    imshow(diff_DoLP_1, []);
    colormap(jet); colorbar;
    caxis([0 1]);
    title(sprintf('|%s - %s|', method1, method2));
    
    subplot(2,2,3);
    imshow(diff_DoLP_2, []);
    colormap(jet); colorbar;
    caxis([0 1]);
    title(sprintf('|%s - %s|', method1, method3));
    
    subplot(2,2,4);
    imshow(diff_DoLP_fourier, []);
    colormap(jet); colorbar;
    caxis([0 1]);
    title(sprintf('|%s - %s|', method2, method3));
end