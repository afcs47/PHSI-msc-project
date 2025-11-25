%% HSI Polarization Results Comparison Only
close all; clear; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi')); %for dolp and aolp comparison functions

%% Select folder containing saved results
resultsFolder = uigetdir(pwd, 'Select Folder Containing Saved Result .mat Files');
if resultsFolder == 0
    disp('No folder selected. Exiting.');
    return
end

%% Output Folder
outputFolder = uigetdir(pwd, 'Select Output Folder for Figures');
if outputFolder == 0
    disp('No output folder selected. Exiting.');
    return
end

[~, folderName] = fileparts(resultsFolder);
outputFolder = fullfile(outputFolder, [folderName ' Methods comparison']); % avoid scattered images
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end % Create the subfolder based on the filename variable if it doesn't exist

%% Extract selectedType + selectedDay from folder name
parts = split(folderName, '_');

selectedType = parts{1};
if numel(parts) >= 2
    selectedDay = parts{2};
else
    selectedDay = "";
end

fprintf('Detected sample type: %s | day: %s\n', selectedType, selectedDay);

%% Load all .mat analysis files

files = dir(fullfile(resultsFolder, '*.mat'));
if isempty(files)
    error('No .mat files found in %s.', resultsFolder);
end

% Load each file into a struct with the same base name
analysisData = struct();
for k = 1:length(files)
    filePath = fullfile(resultsFolder, files(k).name);
    [~, baseName, ~] = fileparts(files(k).name);

    % Load file
    dataStruct = load(filePath);

    % Store inside main struct using filename as field
    %analysisData.(baseName) = dataStruct;

    % Remove any trailing like '_20251122001430'
    cleanedName = regexprep(baseName, '(_\d{8}_\d{6})$', '');

    % Store using cleaned name
    analysisData.(cleanedName) = dataStruct;
end

datasetNames = fieldnames(analysisData);


%% Method comparison across all available methods

disp('Starting method comparison across all polarization methods...');

methodStructs = struct();

% Collect all method results into unified structure
for i = 1:numel(datasetNames)

    fname = datasetNames{i};
    data = analysisData.(fname);

    % Identify method
    method = identifyPolMethodUsed(fname);

    % Store only polarization-relevant fields
    entry = struct();

    if isfield(data, 'S0'),entry.S0 = data.S0; end
    if isfield(data, 'S1'),entry.S1 = data.S1; end
    if isfield(data, 'S2'),entry.S2 = data.S2; end
    if isfield(data, 'S3'),entry.S3 = data.S3; end
    if isfield(data, 'DoLP'), entry.DoLP = data.DoLP; end
    if isfield(data, 'AoLP'),entry.AoLP = data.AoLP; end

    % Save only methods with valid polarization data
    if isfield(entry, 'DoLP')
        methodStructs.(method) = entry;
    end
end


%% Full comparison of polarization paremeters per method
disp('Running detailed polarization-parameter comparisons per method...');

for k = 1:numel(methodNames)
    m = methodNames{k};
    entry = methodStructs.(m);

    s0 = entry.S0;  s1 = entry.S1;  s2 = entry.S2;
    if isfield(entry,'S3'), s3 = entry.S3; else, s3 = []; end
    dolp = entry.DoLP;
    aolp = entry.AoLP;

    plot_pol_parameters_comparison(s0, s1, s2, dolp, aolp, selectedType, m, outputFolder);
end


%% 3-method at a time comparison 

% 1. Standard vs LSQ vs FFT
plotComparisonFor3Methods(methodStructs, 'Standard', 'Fourier_LSQ', 'Fourier_FFT', outputFolder, sprintf('%s_%s_Std_LSQ_FFT_Comparison', selectedType, selectedDay));

% 2. Standard vs LSQ vs Full LSQ
plotComparisonFor3Methods(methodStructs, 'Standard', 'Fourier_LSQ', 'Fourier_Full_LSQ', outputFolder, sprintf('%s_%s_Std_LSQ_FullLSQ_Comparison', selectedType, selectedDay));

% 3. Standard vs FFT vs Full FFT
plotComparisonFor3Methods(methodStructs, 'Standard', 'Fourier_FFT', 'Fourier_Full_FFT', outputFolder,sprintf('%s_%s_Std_FFT_FullFFT_Comparison', selectedType, selectedDay));

% 4. Standard vs SPIE vs SPIE Simple
plotComparisonFor3Methods(methodStructs, 'Standard', 'SPIE', 'SPIE_Simple', outputFolder, sprintf('%s_%s_Std_SPIE_SPIESimple_Comparison', selectedType, selectedDay));

% 5. Standard vs LSQ vs SPIE Simple
plotComparisonFor3Methods(methodStructs, 'Standard', 'Fourier_LSQ', 'SPIE_Simple', outputFolder, sprintf('%s_%s_Std_LSQ_SPIESimple_Comparison', selectedType, selectedDay));


disp('All method comparisons completed successfully.');
