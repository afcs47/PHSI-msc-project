function export_figure(figHandle, filenameBase, outputDir)
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    % Construct filename with timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filenameBase = sprintf('%s_%s.mat', filenameBase, timestamp);

    exportgraphics(figHandle, fullfile(outputDir, filenameBase + ".png"), 'Resolution', 400);
    %exportgraphics(figHandle, fullfile(outputDir, filenameBase + ".eps"), 'ContentType', 'vector');
    %saveas(figHandle, fullfile(outputDir, filenameBase + ".fig"));
end