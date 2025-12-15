
function save_analysis_data(dataStruct, dataType, datasetname, outputDir) 
    % Ensure output directory exists
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Construct filename with timestamp
    timestamp = string(datetime("now","Format","yyyyMMdd_HHmmss"));
    filename = sprintf('%s_%s_%s.mat', dataType, datasetname, timestamp);
    save(fullfile(outputDir, filename), '-struct', 'dataStruct');

    fprintf('Saved %s data for "%s" to %s\n', dataType, datasetname, filename);
end
