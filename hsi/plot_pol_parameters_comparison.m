function plot_pol_parameters_comparison(S0, S1, S2, DoLP, AoP_deg, datasetname, method, outputFolder)
    % Plot Results: S0, S1, S2, DoLP and AoLP (grayscale)

    titleText = [method, ' Polarization Parameters Comparison for ', datasetname];
    % Stretch between 1st and 99th percentiles to brighten the image
    fig = figure('Name', titleText); 
    subplot(3,2,1);
    imshow(uint8(255 * mat2gray(S0, [prctile(S0(:), 1), prctile(S0(:), 99)])), []); title('S0 (Total Intensity)');
    subplot(3,2,2);
    imshow(uint8(255 * mat2gray(S1, [prctile(S1(:), 1), prctile(S1(:), 99)])), []); title('S1');
    subplot(3,2,3);
    imshow(uint8(255 * mat2gray(S2, [prctile(S2(:), 1), prctile(S2(:), 99)])), []); title('S2');
    
    % DoLP grayscale
    subplot(3,2,4);
    imshow(uint8(255 * mat2gray(DoLP, [prctile(DoLP(:), 1), prctile(DoLP(:), 99)]))); title('DoLP');
    
    % AoP grayscale
    subplot(3,2,5);
    imshow(uint8(255 * mat2gray(AoP_deg))); title('BW AoP (Degrees)');
    
    %Smoothing
    kernel = fspecial('average', [2 2]); 
    bw_AoLP_smooth = imfilter(uint8(255 * mat2gray(AoP_deg)), kernel);
    
    subplot(3,2,6);
    imshow(bw_AoLP_smooth); title('Smoothed BW AoP (Degrees)');

    export_figure(fig, [strrep(titleText, ' ', '_')], outputFolder);

end
