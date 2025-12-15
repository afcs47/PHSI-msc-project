%% Pol analysis of HSI files - NOVEMBER 2025

% How the code is structured:
% - Lets the user choose which dataset to load (stokesData for the standard method or fourierData for Fourier methods)
% - Lets the user pick which polarization computation method to run (Fourier LSQ, Fourier Full LSQ, Fourier FFT, Fourier Full FFT, SPIE, SPIE simplified)
% - Processes one method at a time, saving results and clearing memory in between to prevent crashes (Matlab's "Out of memory" error)
% - Always computes reflectances first (needed for phsi setup 3/ fusion approach), then DoLP and AoLP depending on the chosen method
% - Ensures DoLP and AoLP maps are 2D, not 1D!

%Note: Comment lines 159 and 164 from main script, and 29 and 34 from run_standard_polarization_analysis to reduce computation time by not saving all spectral reflectances

%% Some initial settings

close all   % Closes all open figure windows
clear % Clears variables from the workspace
clc % Clears the Command Window

plotComparison = true;

% Get auxiliary matlab functions previously created for HSI data analysis
addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi'));

%% Select parent directory (where datasets are stored)
%baseDir = uigetdir('', 'Select Parent Directory Containing Data Folders');
baseDir = "D:\afili\Transferências\Tese - files\*";
basePath = fileparts(baseDir); % Get path from directory

outputFolderOG = uigetdir(pwd, 'Select Output Folder'); % Ask where to save results
if outputFolderOG == 0
    disp('No output folder selected. Exiting.');
    return
end

%% Load all HSI data
%processHsTestData(baseDir, true); % Confirm dataset separation by displaying the parsed groups
process_hs_testData(baseDir);

%% User decides which type of dataset to load (Standard or Fourier or no WG)
datasetChoice = questdlg('Select data type to process:', 'Dataset Selection', 'Standard (WG 0°, 45°, 90°, 135°)', 'Fourier (WG 0–180°)', 'No polarizer (woWGpolData)', 'Fourier (WG 0–180°)'); % standard for a quicker analysis, fourier for other methods and no polarizer for just reflectances computation

switch datasetChoice
    case 'Standard (WG 0°, 45°, 90°, 135°)'
        load('allHsTestData.mat', 'stokesData'); % Load hsi dataset list (add ',true' to confirm dataset separation by displaying the parsed groups)
        if ~exist('stokesData', 'var'), error('stokesData not found in allTestData.mat'); end % Check data loading
        dataTable = stokesData;

        disp('Running standard 4-angle polarization analysis...');
        run_standard_polarization_analysis(dataTable, basePath, outputFolderOG);
        %return; % Exit after this analysis (no Fourier methods needed)

        %Do NOT exit the script (restart file at output file selection)
        disp('Standard method finished.');
        userChoice = questdlg('Standard method done. Select next action:', 'Continue', 'Choose New Method/Dataset', 'Exit', 'Choose New Method/Dataset');

        if strcmp(userChoice, 'Exit')
            disp('Exiting.');
            return;
        else
            % Restart script cleanly
            clear; close all; clc;
            run(mfilename);
            return;
        end

    case 'Fourier (WG 0–180°)'
        load('allHsTestData.mat', 'fourierData'); % Load hsi dataset list (add ',true' to confirm dataset separation by displaying the parsed groups)
        if ~exist('fourierData', 'var'), error('fourierData not found in allHsTestData.mat'); end % Check data loading
        dataTable = fourierData;

    case 'No polarizer (woWGpolData)'
        load('allHsTestData.mat', 'woWGpolData'); % Load hsi dataset list (add ',true' to confirm dataset separation by displaying the parsed groups)
        if ~exist('woWGpolData', 'var'), error('woWGpolData not found in allHsTestData.mat'); end % Check data loading        
        dataTable = woWGpolData;
        disp('Running no-polarizer (woWGpolData) analysis...');
        run_woWGpol_analysis(dataTable, basePath, outputFolderOG);
        %return; % Exit after finishing this part

        %Do NOT exit the script (restart at dataset selection)
        disp('woWG analysis finished.');
        userChoice = questdlg('Analysis done. Select next action:', 'Continue', 'Choose New Method/Dataset', 'Exit', 'Choose New Method/Dataset');
        
        if strcmp(userChoice, 'Exit')
            disp('Exiting.');
            return;
        else
            clear; close all; clc;
            run(mfilename);
            return;
        end
end

%% MAIN LOOP for fourier case: USER CHOOSES SAMPLE AND METHOD
while true
    
    %% Select sample type (and day, when relevant)
    [filteredData, selectedType, selectedDay] = filter_data_by_sample_and_day(dataTable);

    % Save reflectances
    outFolder = fullfile(outputFolderOG, [selectedType '_' selectedDay]);
    if ~exist(outFolder, 'dir'); mkdir(outFolder); end % Create the subfolder based on the filename variable if it doesn't exist

    %% Compute Reflectances first

    %% Check if reflectances have been computed and saved in a recent file
    disp('Checking for existing reflectance file...');
    % Pattern for file name that accounts for date/time suffix:
    reflectPattern = fullfile(outFolder, ['Reflectances_' selectedType '_' selectedDay '*.mat']);
    
    reflectFiles = dir(reflectPattern);
    useCachedReflectances = false;
    
    % Check if file exists and its age
    if ~isempty(reflectFiles)
        % Find newest reflectance file
        [~, newestIdx] = max([reflectFiles.datenum]);
        newestFile = reflectFiles(newestIdx);
        reflectFile = fullfile(newestFile.folder, newestFile.name);
    
        % Compute file age
        fileAgeHours = (now - newestFile.datenum) * 24;
    
        if fileAgeHours < 24 %one day to consider data as recent - change 
            disp(['Found reflectance file updated ' num2str(fileAgeHours, '%.2f') ' hours ago.']);
            disp('Loading reflectances instead of recomputing...');
            
            loadedData = load(reflectFile);
            mean_ref_struct = loadedData.Spatial_reflectances;
            wavelengths = loadedData.wavelengths;

            useCachedReflectances = true;
        else
            disp('Reflectance file older than 24h: recomputing.');
        end
    else
        disp('No reflectance file found: computing reflectances.');
    end

    %% If no recent reflectances exist => compute them
    if ~useCachedReflectances
    
        disp('Computing reflectances for all angles...');
        % Sort by angle
        anglesDeg = cellfun(@str2double, filteredData.WGpolAngle);
        [sortedAngles, sortIdx] = sort(anglesDeg);
        filteredData = filteredData(sortIdx, :);
    
        mean_ref_struct = struct();
        spectral_ref_struct = struct();
        for i = 1:numel(sortedAngles)
            fpath = fullfile(basePath, filteredData.FolderName{i}, 'capture\');
            fname = filteredData.FolderName{i};
    
            [Data, White, Dark, wavelengths] = read_data(fpath, fname); % Load raw HSI data
            HS_calib = apply_calibration(Data, White, Dark); % Apply radiometric calibration
            hsfiltered = apply_sg_filter(HS_calib); % Apply spectral smoothing (Savitzky–Golay filter)
            mean_ref = mean(hsfiltered, 3); % Compute mean reflectance over wavelengths
            
            mean_ref_struct.(sprintf('angle_%d', sortedAngles(i))) = mean_ref; % Save result by angle
            spectral_ref_struct.(sprintf('angle_%d', sortedAngles(i))) = single(hsfiltered); % reduce size
        end
    
        % Resize all reflectance images to match dimensions - avoid the problem of the angle arrays having different/incompatible sizes between them
        mean_ref_struct = resize_reflectances(mean_ref_struct);
        cube = resize_hs_cubes(spectral_ref_struct); % also average across angles to reduce memory occupied
    
        save_analysis_data(struct('Spectral_reflectances', cube,'Spatial_reflectances', mean_ref_struct, 'wavelengths', wavelengths), 'Reflectances', [selectedType '_' selectedDay], outFolder);
    
        disp('Reflectances computed and saved.');
    end

    %%  Build angle-intensity cube (I_theta)
    
    disp('Building angle-intensity cube (I_theta)...');
    
    fields = fieldnames(mean_ref_struct);
    numAngles = numel(fields);
    [rows, cols] = size(mean_ref_struct.(fields{1}));
    I_theta = zeros(rows, cols, numAngles);
    sortedAngles = zeros(numel(fields),1);
    for i = 1:numAngles
        I_theta(:,:,i) = mean_ref_struct.(fields{i});
        parts = split(fields{i}, '_');
        sortedAngles(i) = str2double(parts{2}); %rebuild sortedAngles variable
    end
    
    disp('I_theta successfully constructed.');

    %% Choose Fourier method

    methodOptions = { ...
    'Fourier LSQ (S0,S1,S2)', ...
    'Fourier Full LSQ (S0,S1,S2,S3)', ...
    'Fourier FFT (S0,S1,S2)', ...
    'Fourier Full FFT (S0,S1,S2,S3)', ...
    'SPIE-based Fourier', ...
    'Simplified SPIE-based Fourier'};

    [methodIdx, tf] = listdlg('PromptString', 'Select polarization method to compute:', 'SelectionMode', 'single', 'ListString', methodOptions, 'ListSize', [300, 220]); % Wider dialog for readability


    if ~tf
    disp('No method selected. Exiting.');
    return;
    end

    methodChoice = methodOptions{methodIdx};

    %% Process data by selected method

    switch methodChoice
        %% FOURIER LSQ (Linear)
        case 'Fourier LSQ (S0,S1,S2)'
            [S0, S1, S2, DoLP, AoLP] = compute_fourier_pol(I_theta, sortedAngles, numAngles);
            [S0, S1, S2, DoLP, AoLP] = reshape_variables(S0, S1, S2, DoLP, AoLP, rows, cols);
            
            plot_polarization_spatially(DoLP, 'jet', ['DoLP Fourier LSQ - ' selectedType], outFolder);
            plot_polarization_spatially(AoLP, 'hsv', ['AoLP Fourier LSQ - ' selectedType], outFolder);

            save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'DoLP', DoLP, 'AoLP', AoLP), 'Fourier_LSQ_Results', selectedType, outFolder);

        %% FOURIER FULL LSQ (in case of QWP present in the setup)
        case 'Fourier Full LSQ (S0,S1,S2,S3)'
            [S0, S1, S2, S3, DoLP, AoLP, Ellipt] = compute_fourier_pol_full(I_theta, sortedAngles, numAngles);
            [S0, S1, S2, DoLP, AoLP] = reshape_variables(S0, S1, S2, DoLP, AoLP, rows, cols);

            plot_polarization_spatially(DoLP, 'jet', ['DoLP Full LSQ - ' selectedType], outFolder);
            plot_polarization_spatially(AoLP, 'hsv', ['AoLP Full LSQ - ' selectedType], outFolder);

            save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'S3', S3, 'DoLP', DoLP, 'AoLP', AoLP, 'Ellipticity', Ellipt), 'Fourier_Full_LSQ_Results', selectedType, outFolder);

        %% FOURIER FFT (Linear)
        case 'Fourier FFT (S0,S1,S2)'
            [S0, S1, S2, DoLP, AoLP] = compute_fourier_pol_fft_linear(I_theta, numAngles);
            [S0, S1, S2, DoLP, AoLP] = reshape_variables(S0, S1, S2, DoLP, AoLP, rows, cols);

            plot_polarization_spatially(DoLP, 'jet', ['DoLP FFT - ' selectedType], outFolder);
            plot_polarization_spatially(AoLP, 'hsv', ['AoLP FFT - ' selectedType], outFolder);

            save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'DoLP', DoLP, 'AoLP', AoLP), 'Fourier_FFT_Results', selectedType, outFolder);

        %% FOURIER FULL FFT (in case of QWP present in the setup)
        case 'Fourier Full FFT (S0,S1,S2,S3)'
            [S0, S1, S2, S3, DoLP, AoLP, Ellipt] = compute_fourier_pol_full_fft(I_theta, numAngles);
            [S0, S1, S2, DoLP, AoLP] = reshape_variables(S0, S1, S2, DoLP, AoLP, rows, cols);

            plot_polarization_spatially(DoLP, 'jet', ['DoLP Full FFT - ' selectedType], outFolder);
            plot_polarization_spatially(AoLP, 'hsv', ['AoLP Full FFT - ' selectedType], outFolder);

            save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'S3', S3, 'DoLP', DoLP, 'AoLP', AoLP, 'Ellipticity', Ellipt), 'Fourier_Full_FFT_Results', selectedType, outFolder);

        %% SPIE-BASED 
        case 'SPIE-based Fourier'
            [S0, S1, S2, DoLP, AoLP] = compute_spie(I_theta);
            %[S0, S1, S2, DoLP, ~] = reshape_variables(S0, S1, S2, DoLP,AoLP, rows, cols); % AoLP already in degrees

            plot_polarization_spatially(DoLP, 'jet', ['DoLP SPIE - ' selectedType], outFolder);
            plot_polarization_spatially(AoLP, 'hsv', ['AoLP SPIE - ' selectedType], outFolder);

            save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'DoLP', DoLP, 'AoLP', AoLP), 'SPIE_based_Results', selectedType, outFolder);

        %% SIMPLIFIED SPIE-BASED 
        case 'Simplified SPIE-based Fourier'
            [S0, S1, S2, DoLP, AoLP] = compute_spie_simplified(I_theta);
            %[S0, S1, S2, DoLP, ~] = reshape_variables(S0, S1, S2, DoLP, AoLP, rows, cols); % AoLP already in degrees

            plot_polarization_spatially(DoLP, 'jet', ['DoLP SPIE Simple - ' selectedType], outFolder);
            plot_polarization_spatially(AoLP, 'hsv', ['AoLP SPIE Simple - ' selectedType], outFolder);

            save_analysis_data(struct('S0', S0, 'S1', S1, 'S2', S2, 'DoLP', DoLP, 'AoLP', AoLP), 'SPIE_Simplified_Results', selectedType, outFolder);
    end

    %% Optional comparison
    if plotComparison
        disp('Plotting comparisons between computed maps...');
        plot_pol_parameters_comparison(S0, S1, S2, DoLP, AoLP, selectedType, methodChoice, outFolder);
    end

    %% Clean memory before next method to avoid "Out of memory" Matlab error
    disp('Clearing large variables to free memory...');
    clear Data White Dark HS_calib hsfiltered DoLP AoLP S0 S1 S2 S3 Ellipt I_theta mean_ref_struct;
    
    choice = questdlg('Run another dataset for the same method (Fourier (WG 0–180°))?', 'Next Step', 'Yes', 'Change method', 'Exit', 'Exit');
    if strcmp(choice, 'Exit')
        disp('Processing complete. Exiting.');
        break;
    elseif strcmp(choice, 'Change method')
        disp('Restarting selection...');
        clear dataTable; close all; clc;
        run(mfilename); % restart script
    end
end

