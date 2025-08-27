%% Plot saved HSI + Polarization analysis results
close all;
clear;
clc;

addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi')); % Functions folder
addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi'));


% Ask for folder with saved results
resultsFolder = uigetdir(pwd, 'Select Folder Containing Saved Analysis Data');
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

%% Load saved analysis data
% Expecting files created by save_analysis_data():
%   Standard_Stokes_Results_*.mat
%   Fourier_Stokes_Results_*.mat
%   SPIE_Fourier_Stokes_Results_*.mat
%   SpieSimple_Fourier_Stokes_Results_*.mat

files = dir(fullfile(resultsFolder, '*.mat'));

% Load each file into a struct with the same base name
for k = 1:length(files)
    [~, name] = fileparts(files(k).name);
    dataStruct = load(fullfile(resultsFolder, files(k).name));
    analysisData.(name) = dataStruct;
end

%% Plot each dataset
sampleNames = fieldnames(analysisData);

for i = 1:numel(sampleNames)
    name = sampleNames{i};
    data = analysisData.(name);
    name = strrep(name, '_', ' ');
    
    fprintf('Plotting results for: %s\n', name);

    % Standard results
    if isfield(data, 'Standard_DoLP_map')
        plot_polarization_spatially(data.Standard_DoLP_map, 'jet', ['Standard DoLP - ' name], outputFolder);
        plot_polarization_spatially(data.Standard_AoLP_map, 'hsv', ['Standard AoLP - ' name], outputFolder);
    end

    % Fourier results
    if isfield(data, 'Fourier_DoLP_map')
        plot_polarization_spatially(data.Fourier_DoLP_map, 'jet', ['Fourier DoLP - ' name], outputFolder);
        plot_polarization_spatially(data.Fourier_AoLP_map, 'hsv', ['Fourier AoLP - ' name], outputFolder);
    end

    % SPIE Fourier results
    if isfield(data, 'Spie_Fourier_DoLP_map')
        plot_polarization_spatially(data.Spie_Fourier_DoLP_map, 'jet', ['SPIE Fourier DoLP - ' name], outputFolder);
        plot_polarization_spatially(data.Spie_Fourier_AoLP_map, 'hsv', ['SPIE Fourier AoLP - ' name], outputFolder);
    end

    % SPIE simplified results
    if isfield(data, 'SpieSimple_Fourier_DoLP_map')
        plot_polarization_spatially(data.SpieSimple_Fourier_DoLP_map, 'jet', ['SPIE Simple DoLP - ' name], outputFolder);
        plot_polarization_spatially(data.SpieSimple_Fourier_AoLP_map, 'hsv', ['SPIE Simple AoLP - ' name], outputFolder);
    end

    % Reflectances
    if isfield(data, 'Spatial_reflectances')
        plot_reflectances_spatially(mean(data.Spatial_reflectances,3), 'gray', ['Spatial Reflectances (mean) - ' name], outputFolder);
        plot_reflectances_spatially(mean(data.Spatial_reflectances,3), 'jet', ['Spatial Reflectances (mean) - ' name], outputFolder);
    end
    if isfield(data, 'Wavelength_reflectances')
        mean_ref = mean(data.Wavelength_reflectances,3);
        plot_reflectances_spatially(mean_ref, 'gray', ['Normalized Reflectance - ' name ' (grayscale)' ], outputFolder);
        plot_reflectances_spatially(mean_ref, 'jet', ['Normalized Reflectance - ' name], outputFolder);

        lambda = 550; % example wavelength
        [~, idx] = min(abs(data.wavelengths - lambda));
        plot_reflectances_spatially(data.Wavelength_reflectances(:,:,idx), 'jet', sprintf('Reflectance @ %.0f nm - %s', lambda, name), outputFolder);
    end
end

disp('All plots generated successfully.');
