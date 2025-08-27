function plot_reflectances_spatially(reflectance_map, cmap, titleStr, outputFolder)
    figure("Name", titleStr);
    imshow(reflectance_map, []);
    axis image;
    colormap(cmap);
    colorbar;
    title(titleStr);
    caxis([0 1]);
    
    % Save figure
    saveas(gcf, fullfile(outputFolder, [regexprep(titleStr,'\s+','_') '.png']));
end