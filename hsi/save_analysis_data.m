
function save_analysis_data(dataStruct, dataType, datasetname, outputDir)
    % Ensure output directory exists
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Construct filename with timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('%s_%s_%s.mat', dataType, datasetname, timestamp);
    save(fullfile(outputDir, filename), '-struct', 'dataStruct');

    fprintf('Saved %s data for "%s" to %s\n', dataType, datasetname, filename);
end
