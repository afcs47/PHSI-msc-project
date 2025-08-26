function processHsTestData(baseDir, verifyDatasetSeparation) % Parses and categorizes hyperspectral data folders into structured datasets based on the given foulder names, optionally displaying parsed groups info

    if nargin < 2
        verifyDatasetSeparation = false;
    end

    % Get subfolders
    folderList = dir(baseDir);
    folderList = folderList([folderList.isdir]);
    folderList = folderList(~ismember({folderList.name}, {'.', '..'}));

    % Initialize result array
    results = {};

    % Pattern for test folders
    testPattern = 'test_([a-zA-Z]+)(.+?)_(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})';

    % Parse folders
    for i = 1:length(folderList)
        fname = folderList(i).name;

        tokens = regexp(fname, testPattern, 'tokens');

        row = {fname, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A'}; % Default row

        if ~isempty(tokens)
            tokens = tokens{1};
            row{2} = tokens{1}; % SampleName
            row{3} = tokens{2}; % SampleDetails
            row{4} = tokens{3}; % Date
            row{5} = tokens{4}; % Time

            % Extract day from SampleDetails
            dayMatch = regexp(tokens{2}, '(\d+)[wW]?[gG]?', 'tokens');
            if ~isempty(dayMatch)
                row{6} = dayMatch{1}{1};
            end

            % Extract angle
            angleMatch = regexp(tokens{2}, '(\d+)(deg)', 'tokens');
            if ~isempty(angleMatch)
                row{7} = angleMatch{1}{1};
            end
        end

        % Add 'pol' folder manually
        if startsWith(fname, 'pol')
            row = {fname, 'N/A', 'jun polarization', 'N/A', 'N/A', 'N/A', 'N/A'};
        end

        results = [results; row];
    end

    % Convert to table
    resultsTable = cell2table(results, ...
        'VariableNames', {'FolderName', 'SampleName', 'SampleDetails', 'Date', 'Time', 'Day', 'WGpolAngle'});

    % Keep only most recent files by SampleName + SampleDetails
    [~, ia] = unique(resultsTable(:, 2:3), 'last');
    resultsTable = resultsTable(ia, :);
    
    if verifyDatasetSeparation
        %size(resultsTable)
        disp(resultsTable);
    end

    %% Dataset Grouping
    testRows = startsWith(resultsTable.FolderName, 'test_');
    polRows  = startsWith(resultsTable.FolderName, 'pol');

    testData = resultsTable(testRows, :);
    polData  = resultsTable(polRows, :);

    % Associate pol folders to test data based on Day
    testData.PolFolder = repmat({'N/A'}, height(testData), 1); % preallocate column
    for i = 1:height(testData)
        day = testData.Day{i};
        polMatch = find(contains(polData.FolderName, [day 'jun']), 1);
        if ~isempty(polMatch)
            testData.PolFolder{i} = polData.FolderName{polMatch};
        end
    end

    % Categorize
    stokesData = testData(ismember(testData.WGpolAngle, {'0','45','90','135'}), :);
    fourierData = testData(~ismember(testData.WGpolAngle, {'N/A','woWGpol','woWGpolcalib'}), :);
    woWGpolData = testData(contains(testData.SampleDetails, 'woWGpol') & ~contains(lower(testData.SampleDetails), 'calib'), :);
    calibHsData = testData(strcmp(testData.WGpolAngle, 'N/A') & contains(lower(testData.SampleDetails), 'calib'), :);

    % Groupings
    groupedStokes = groupsummary(stokesData, 'SampleDetails');
    groupedFourier = groupsummary(fourierData, 'SampleDetails');
    groupedSamples = groupsummary(testData, 'SampleName');

    if verifyDatasetSeparation
        disp('--- Stokes Data ---'); disp(stokesData);
        disp('--- Stokes Groups ---'); disp(groupedStokes);
        disp('--- Fourier Data ---'); disp(fourierData);
        disp('--- Fourier Groups ---'); disp(groupedFourier);
        disp('--- Calibration Files ---'); disp(calibHsData);
        disp('--- woWGpol Data ---'); disp(woWGpolData);
        disp('--- Sample Groups ---'); disp(groupedSamples);
    end

    % Save all variables into one .mat file
    save('allHsTestData.mat', 'stokesData', 'fourierData', 'woWGpolData', 'calibHsData', 'groupedStokes', 'groupedFourier', 'groupedSamples', 'testData', 'polData');

    fprintf('Hyperspectral Data successfully processed and saved to allTestData.mat\n');
end
