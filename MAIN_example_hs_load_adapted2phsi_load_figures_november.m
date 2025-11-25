%% Plot saved HSI + Polarization analysis results
close all;
clear;
clc;

addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi')); % Functions folder
%addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi'));


% Ask for folder with saved results
resultsFolder = uigetdir(pwd, 'Select Folder Containing the Saved Analysis Data');
if resultsFolder == 0
    disp('No folder selected. Exiting.');
    return
end

% Ask output folder for saving plots
outputFolder = uigetdir(pwd, 'Select Output Folder for Figures');
if outputFolder == 0
    disp('No output folder selected. Exiting.');
    return
end

[~, folderName] = fileparts(resultsFolder);
outputFolder = fullfile(outputFolder, [folderName ' Reloaded plots']); % avoid scattered images
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end % Create the subfolder based on the filename variable if it doesn't exist


%% Load saved analysis data
files = dir(fullfile(resultsFolder, '*.mat'));
if isempty(files)
    error('No .mat files found in selected results folder.');
end

%% Extract the part between last underscore and ".mat"

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

%% Plot each dataset
datasetNames = fieldnames(analysisData);

% Loop through all loaded datasets
for i = 1:numel(datasetNames)
    fname = datasetNames{i};
    data = analysisData.(fname);
    fprintf('Plotting results for: %s\n', fname);
    fname = strrep(fname, '_', ' ');

    %% Spatial Reflectances
    if isfield(data, 'Spatial_reflectances') % Plot mean reflectance

        spatial = data.Spatial_reflectances;

        % Mean reflectance map
        fields = fieldnames(spatial);
        firstImg = spatial.(fields{1});
        meanImg = zeros(size(firstImg));

        for f = 1:numel(fields)
            meanImg = meanImg + spatial.(fields{f});
        end
        meanImg = meanImg ./ numel(fields); %normalization

        plot_reflectances_spatially(meanImg, 'gray', ['Spatial Reflectance (mean) - ' fname], outputFolder);

        plot_reflectances_spatially(meanImg, 'jet', ['Spatial Reflectance (mean) - ' fname], outputFolder);

    end

    %% Polarization Parameters
    
    if isfield(data, 'DoLP')
        plot_polarization_spatially(data.DoLP, 'jet', ['DoLP - ' fname], outputFolder);
    end

    if isfield(data, 'AoLP')
        plot_polarization_spatially(data.AoLP, 'hsv', ['AoLP - ' fname], outputFolder);
    end

    %% Spectral reflectances
    if isfield(data, 'wavelengths') && isfield(data, 'Spectral_reflectances')

        wavelengths = data.wavelengths;

        lambdaList = [400, 450, 550, 650, 700, 780, 1000]; % - Visible light: from ~380nm to ~750nm; Violet~400nm; Blue~450nm; Green~550nm; Red~700nm; - Near-Infrared (IR-A) ~780nm to ~1.4 um;    Effective range of HSI camera: 396.7-1004.2 nm


        fields = fieldnames(data.Spectral_reflectances);
        cube = data.Spectral_reflectances.(fields{1}); % pick first one

        % Plot reflectance at specific wavelengths 
        for lambda = lambdaList
            % Find nearest wavelength index
            [~, idx] = min(abs(wavelengths - lambda));

            img = cube(:,:,idx);

            plot_reflectances_spatially(img, 'jet', sprintf('Normalized Reflectance @ %.0f nm - %s', lambda, fname), outputFolder);
            plot_reflectances_spatially(img, 'gray', sprintf('Normalized Reflectance @ %.0f nm (grayscale) - %s', lambda, fname), outputFolder);
        end
    end

end

disp('All plots generated successfully.');
