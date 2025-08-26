function processPolTestData(basePath, selectedSampleType, verifyMatches) % Parses and categorizes polarization data files based on the given file names, optionally displaying matched file info

    % If verifyMatches is not passed in, default to false
    if nargin < 4
        verifyMatches = false;
    end

    % Build full path to the "pol data" folder
    polFolder = fullfile(basePath, 'pol data');
    %polFolder = fullfile(basePath); % use if basePath already includes pol data directory

    % Check if the folder exists; throw error if not found
    if ~isfolder(polFolder)
        error('"pol data" folder not found at: %s', basePath);
    end

    % List all .raw files in the "pol data" directory
    rawFiles = dir(fullfile(polFolder, '*.raw'));

    % Initialize a cell array to store parsed file info
    results = {};

    % Regular expression pattern to match:
    pattern = '^([a-zA-Z0-9]+)?(_woWG)?(_calib)?_(\d{1,2})([a-zA-Z]{3})\.raw$'; % <sample_type><version>[_woWG][_calib]_<day>.raw

    % Loop over all found .raw files to select the needed files
    for i = 1:length(rawFiles)
        fname = rawFiles(i).name;

        % Apply regex to extract tokens from filename
        tokens = regexp(fname, pattern, 'tokens');

        % Skip if the filename does not match the expected format
        if isempty(tokens)
            continue;
        end

        % Extract matched tokens from regex
        tokens = tokens{1};

        % Assign token values to descriptive variables
        fullSample = tokens{1}; % e.g. 'solutions2'
        isWoWG = ~isempty(tokens{2}); % true if _woWG is present
        isCalib = ~isempty(tokens{3}); % true if _calib is present
        fileDay = tokens{4}; % e.g. '17'

        % Extract sampleType and sampleVersion from fullSample
        sampleMatch = regexp(fullSample, '^([a-zA-Z]+)(\d*)$', 'tokens');

        if ~isempty(sampleMatch)
            sampleMatch = sampleMatch{1}; % 1×1 cell array '{1×2 cell}' -> 1×2 cell array '{'solutions'}    {'2'}'
            sampleType = sampleMatch{1}; % 'solutions'
            sampleVersion = sampleMatch{2};  % '2' (can be empty)
        else
            sampleType = fullSample; 
            sampleVersion = '';
        end
    
        % Fill in defaults if needed
        if isempty(sampleVersion)
            sampleVersion = '1';
        end


        % Determine category of file based on flags
        if isWoWG && isCalib
            category = "woWGCalib";
        elseif isWoWG
            category = "woWG"; 
        elseif isCalib
            category = "Calib"; 
        else
            category = "Standard"; 
        end

        % Add parsed data to results array
        row = { ...
            fname, ...  % FileName
            sampleType, ... 
            sampleVersion, ... 
            isWoWG, ...  %(true/false)
            isCalib, ...  %(true/false)
            fileDay, ...  
            category ... %(Standard/Calib/woWG/woWGCalib)
        };

        % Append row to results
        results = [results; row];
    end

    % Convert cell array to table for structured handling
    resultsTable = cell2table(results, 'VariableNames', {'FileName', 'SampleType', 'SampleVersion', 'IsWoWG', 'IsCalib', 'Day', 'Category'});

    % Filter table by selected sample and day
    matchIdx = strcmpi(resultsTable.SampleType, selectedSampleType);
    matchedFiles = resultsTable(matchIdx, :);  % Only keep matching files

    % Split matched files into categorized tables
    standardData = matchedFiles(strcmp(matchedFiles.Category, 'Standard'), :);
    calibPolData = matchedFiles(strcmp(matchedFiles.Category, 'Calib'), :);
    woWGData = matchedFiles(strcmp(matchedFiles.Category, 'woWG'), :);
    woWGCalibData = matchedFiles(strcmp(matchedFiles.Category, 'woWGCalib'), :);

    % If requested, display all filtered results and groups
    if verifyMatches
        disp('--- Matched Files ---'); disp(matchedFiles);
        disp('--- Standard Files ---'); disp(standardData);
        disp('--- Calib Files ---'); disp(calibPolData);
        disp('--- woWG Files ---'); disp(woWGData);
        disp('--- woWG + Calib Files ---'); disp(woWGCalibData);
    end

    % Save all categorized results to a .mat file for future use
    save('processedPolRawFiles.mat', 'matchedFiles', 'standardData', 'calibPolData', 'woWGData', 'woWGCalibData');

    % Print final summary to console
    fprintf('Processed %d .raw files for "%s".\n', height(matchedFiles), selectedSampleType);

end
