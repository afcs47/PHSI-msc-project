function save_polarization_results(outputFolder, base_name,  DoLP, AoLP, DoLP_2nd, AoLP_2nd, DoLP_4th, AoLP_4th, bw_DoLP, bw_AoLP, bw_DoLP_smooth, bw_AoLP_smooth)

    % Ensure output path exists
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Save MAT data for later analysis
    save(fullfile(outputFolder, [base_name '_polarization_data.mat']),  'DoLP', 'AoLP', 'DoLP_2nd', 'AoLP_2nd', 'DoLP_4th', 'AoLP_4th', 'bw_DoLP', 'bw_AoLP', 'bw_DoLP_smooth', 'bw_AoLP_smooth');

    % Save grayscale images
    imwrite(bw_DoLP, fullfile(outputFolder, [base_name '_DoLP_2nd_bw.png']));
    imwrite(bw_AoLP, fullfile(outputFolder, [base_name '_AoLP_2nd_bw.png']));
    imwrite(bw_DoLP_smooth, fullfile(outputFolder, [base_name '_DoLP_2nd_bw_smooth.png']));
    imwrite(bw_AoLP_smooth, fullfile(outputFolder, [base_name '_AoLP_2nd_bw_smooth.png']));

    % Optionally, save AoLP in color (HSV-mapped RGB)
    rgb_AoLP = map_aolp_to_rgb(mod(mean(AoLP_2nd,3), pi), 0, 180);
    imwrite(rgb_AoLP, fullfile(outputFolder, [base_name '_AoLP_2nd_rgb.png']));

    fprintf('All results saved to: %s\n', outputFolder);
end
