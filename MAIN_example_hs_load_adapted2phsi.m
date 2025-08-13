
%hsi + pol (multiple datasets and samples) - just spatial dolp and aolp computation

close all   % Closes all open figure windows
clear % Clears variables from the workspace
clc % Clears the Command Window

plotComparison = true;

% % Get path of current script
% currentScriptPath = fileparts(mfilename('fullpath'));
% % Folder with custom functions
% hsiFuncFolder = fullfile(currentScriptPath, 'hsi');
% % Add it to MATLAB path
% addpath(hsiFuncFolder);

addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi')); % Get auxiliary matlab functions previously created for HSI data analysis

% Select parent directory (where datasets are stored)
%baseDir = uigetdir('', 'Select Parent Directory Containing Data Folders');
baseDir = "D:\afili\Transferências\Tese - files\*";
basePath = fileparts(baseDir); % Get path from directory

outputFolder = uigetdir(pwd, 'Select Output Folder'); % Ask where to save results
if outputFolder == 0
    disp('No output folder selected. Exiting.');
    return
end

%% LOad all HSI data
%processHsTestData(baseDir, true); % Confirm dataset separation by displaying the parsed groups
processHsTestData(baseDir);

% Load data table with valid folders' name
load('allHsTestData.mat', 'stokesData'); % Only loads folders necessary for Stokes parameters' standard calculation (WG at 0, 45, 90 and 135 deg)

% Check data loading
if ~exist('stokesData', 'var')
    error('stokesData not found in allTestData.mat');
end

while true
%% Select sample type
sampleTypes = unique(stokesData.SampleName);
[selectedTypeIndex, tf] = listdlg('PromptString','Select sample type:', 'SelectionMode','single', 'ListString', sampleTypes);
if ~tf, return; end
selectedType = sampleTypes{selectedTypeIndex};

% Filter stokesData by selected type by comparing strings/names
filteredData = stokesData(strcmp(stokesData.SampleName, selectedType), :);

dayList = unique(stokesData.Day(strcmp(stokesData.SampleName, selectedType)));
if length(dayList)>1
    [selectedDayIndex, tf] = listdlg('PromptString','Select sample acquisition day:', 'SelectionMode','single', 'ListString', dayList);
    if ~tf, return; end
    selectedDay = dayList{selectedDayIndex};
    filteredData = filteredData(strcmp(filteredData.Day, selectedDay), :);

else
    selectedDay = dayList{1};
end


%% Create a subfolder based on the filename variable to save the data
datasetFolder = fullfile(outputFolder, [selectedType '_' selectedDay]);
% Create the folder if it doesn't exist
if ~exist(datasetFolder, 'dir')
    mkdir(datasetFolder);
end
% Use datasetFolder as the new output directory
outputFolder = datasetFolder;

%% Generate options from datasets available
options = {};
for i = 1:height(filteredData)
    label = sprintf('%s (%s jun) - %s', selectedType, selectedDay, filteredData.WGpolAngle{i});
    folderName = filteredData.FolderName{i};
    fullPath = fullfile(basePath, folderName, 'capture\');
    angle = filteredData.WGpolAngle{i};
    options(end+1, :) = {label, folderName, fullPath angle};
end

%% Load data acquired without the WG polarizer (HSI only data) corresponding to the selected sample type 

% % Load data table with valid folders' name
% load('allHsTestData.mat', 'woWGpolData'); % Only loads relevant folders
% 
% idx = strcmp(woWGpolData.SampleName, selectedType); %Find HSI only data folder corresponding to the selected sample type
% if sum(idx) ~= 1
%     error('Error: Expected exactly one woWGpolData match for selected type, but got %d.', sum(idx));
% end
% folderName = woWGpolData.FolderName{idx};  % Extract name as string
% woWGPath = fullfile(basePath, folderName);
% 
% if ~isfolder(woWGPath) % Check if folder exists
%     warning('WG-free (HSI only) folder not found: %s', woWGPath);
%     woWGPath = '';  % Handle gracefully 
% end
% 
% 
% % Add hsi only data to options
% label = sprintf('%s (%s jun) - %s', selectedType, woWGpolData.Day{1}, 'wo WG pol'); % Assuming all wo WG pol data was taken on the same day
% fullPath = fullfile(woWGPath, 'capture\');
% angle = 'N/A';
% options(end+1, :) = {label, folderName, fullPath, angle};

% [selectedIndexes, tf] = listdlg('PromptString','Select datasets to process:', 'SelectionMode','multiple', 'ListString', options(:,1));
% if ~tf, return; end


%% Standard computation of DoLP and AoLP
% DoLP and AoLP as a spatial function - Standard method
spatial_reflectances = struct();
for k = 1:length(options)
    datasetname = options{k, 1};
    fname_selected = options{k, 2};
    fpath = options{k, 3};
    angle = options{k, 4};

    % Extract the calibrated hyperspectral data (HS_calibrated) and region of interest (ROI)
    [Data, White, Dark, wavelengths] = read_data(fpath, fname_selected); % Read HS image
    HS_calibrated = apply_calibration(Data, White, Dark); % Calibration
    [R, G, B] = fake_rgb(false, HS_calibrated, datasetname); % Fake RGB image
    hsfiltered = apply_sg_filter(HS_calibrated); % Savitzky-Golay Filtering
    
    % Compute mean reflectance over the wavelengths
    mean_ref = mean(hsfiltered, 3); % x, y, lambda

   
    % Use angle to make valid struct field (needs to start with a letter)
    if strcmp(angle, 'N/A')
        fieldName = 'no_angle';
    else
        fieldName = ['angle_' angle];
    end

    % Store results
    spatial_reflectances.(fieldName) = mean_ref;

end

%% Problem: the angle arrays have different/incompatible sizes between them
spatial_reflectances = resize_reflectances(spatial_reflectances);

%% Compute DoLP for each pixel
DoLP_wg = compute_dolp(spatial_reflectances);
plot_polarization_spatially(DoLP_wg, 'jet', ['Standard Spatial Distribution of DoLP for ' selectedType], outputFolder);

% Compute AoLP for each pixel
AoLP_wg = compute_aolp(spatial_reflectances);
plot_polarization_spatially(AoLP_wg, 'hsv', ['Standard Spatial Distribution of AoLP for ' selectedType], outputFolder);

if plotComparison; plot_pol_parameters_comparison(spatial_reflectances.('angle_0') + spatial_reflectances.('angle_90'), spatial_reflectances.('angle_0') - spatial_reflectances.('angle_90'), spatial_reflectances.('angle_45') - spatial_reflectances.('angle_135'), DoLP_wg, AoLP_wg, selectedType, 'Standard', outputFolder); end

%% Save spatial DoLP/AoLP
standard_data_struct = struct('Standard_DoLP_map', DoLP_wg, 'Standard_AoLP_map', AoLP_wg, 'Stokes_S0', spatial_reflectances.('angle_0') + spatial_reflectances.('angle_90'), 'Stokes_S1', spatial_reflectances.('angle_0') - spatial_reflectances.('angle_90'), 'Stokes_S2', spatial_reflectances.('angle_45') - spatial_reflectances.('angle_135'), 'wavelenghts', wavelengths);

% Save the full struct to file
save_analysis_data(standard_data_struct, 'Standard_Stokes_Results', selectedType, outputFolder);

%% Fourier-based computation of DoLP and AoLP
% Load data table with WG polarizer at 0 - 180 deg in 5 deg steps
load('allHsTestData.mat', 'fourierData');

% Filter data for the selected sample and day
filteredFourier = fourierData(strcmp(fourierData.SampleName, selectedType) & strcmp(fourierData.Day, selectedDay), :);

% Ensure data includes angles from 0 to 180 in steps of 5
anglesDeg = cellfun(@str2double, filteredFourier.WGpolAngle);
[sortedAngles, sortIdx] = sort(anglesDeg); % Sort by angle
filteredFourier = filteredFourier(sortIdx, :); % Apply sorting to table

% Build struct of mean reflectances for each angle
mean_ref_struct = struct();

% Process each angular dataset
for i = 1:numel(sortedAngles)
    fpath = fullfile(basePath, filteredFourier.FolderName{i}, 'capture\');
    fname = filteredFourier.FolderName{i};

    [Data, White, Dark, ~] = read_data(fpath, fname); % Load hyperspectral data and references
    HS_calibrated = apply_calibration(Data, White, Dark); % Perform radiometric calibration
    hsfiltered = apply_sg_filter(HS_calibrated); % Apply Savitzky-Golay filtering
    
    % Compute mean reflectance across all wavelengths
    mean_ref = mean(hsfiltered, 3); % [rows × cols]

    fieldName = sprintf('angle_%d', sortedAngles(i)); % e.g., 'angle_0', 'angle_5', etc.
    mean_ref_struct.(fieldName) = mean_ref;
end

% Resize all reflectance images using existing function
mean_ref_struct = resize_reflectances(mean_ref_struct);

% Extract angles and reshape for Fourier fitting
fields = fieldnames(mean_ref_struct);
numAngles = numel(fields);
theta_vals = zeros(numAngles, 1);
for i = 1:numAngles
    theta_vals(i) = str2double(erase(fields{i}, 'angle_'));
end

% Sort angles and corresponding image stack
[theta_vals, sortIdx] = sort(theta_vals);
[rows, cols] = size(mean_ref_struct.(fields{1}));
I_theta = zeros(rows, cols, numAngles);

for i = 1:numAngles
    I_theta(:, :, i) = mean_ref_struct.(fields{sortIdx(i)});
end

%% Compute Fourier coefficients - simple method
%[~, ~, ~, DoLP_fourier, AoLP_fourier] = compute_fourier_pol(I_theta, sortedAngles, numAngles);
[S0, S1, S2, DoLP_fourier, AoLP_fourier] = compute_fourier_pol(I_theta, sortedAngles, numAngles);

% Reshape to image format
DoLP_fourier_img = reshape(DoLP_fourier, rows, cols);
AoLP_fourier_img = rad2deg(reshape(AoLP_fourier, rows, cols)); % Convert to degrees

S0 = reshape(S0, rows, cols); S1 = reshape(S1, rows, cols); S2 = reshape(S2, rows, cols); 

if plotComparison; plot_pol_parameters_comparison(S0, S1, S2, DoLP_fourier_img, AoLP_fourier_img, selectedType, 'Fourier', outputFolder); end

% Plot Fourier-based DoLP and AoLP
plot_polarization_spatially(DoLP_fourier_img, 'jet', ['Fourier Spatial Distribution of DoLP for ' selectedType], outputFolder);
plot_polarization_spatially(AoLP_fourier_img, 'hsv', ['Fourier Spatial Distribution of AoLP for ' selectedType], outputFolder);

%% Save results
fourier_data_struct = struct('Fourier_DoLP_map', DoLP_fourier_img, 'Fourier_AoLP_map', AoLP_fourier_img, 'Stokes_S0', S0, 'Stokes_S1', S1, 'Stokes_S2', S2, 'wavelenghts', wavelengths);

% Save the full struct to file
save_analysis_data(fourier_data_struct, 'Fourier_Stokes_Results', selectedType, outputFolder);
%% FOURIER ANALYSIS BASED ON SPIE SENSOR MODEL
%[~, ~, ~, DoLP, AoP_deg] = compute_spie(I_theta, true);
[S0, S1, S2, DoLP, AoP_deg] = compute_spie(I_theta);

if plotComparison; plot_pol_parameters_comparison(S0, S1, S2, DoLP, AoP_deg, selectedType, 'SPIE', outputFolder); end

% Display of color-coded DoLP and AoP
plot_polarization_spatially(DoLP, 'jet', ['SPIE Fourier Spatial Distribution of DoLP for ' selectedType], outputFolder);
plot_polarization_spatially(AoP_deg, 'hsv', ['SPIE Fourier Spatial Distribution of AoLP for ' selectedType], outputFolder);

%% Save results
spie_data_struct = struct( 'Spie_Fourier_DoLP_map', DoLP, 'Spie_Fourier_AoLP_map', AoP_deg, 'Stokes_S0', S0, 'Stokes_S1', S1, 'Stokes_S2', S2, 'wavelenghts', wavelengths);

% Save the full struct to file
save_analysis_data(spie_data_struct, 'SPIE_Fourier_Stokes_Results', selectedType, outputFolder);

%% Simplified SPIE paper method (m12 = 0, m13 = 0, Si_int = 0, offset = 0 BUT C1 is also unaccounted for/ unknown) => relative contrast analysis of DoLP and AoLP since Stokes parameters are scaled (unitless)
%[~, ~, ~, DoLP_img, AoP_img] = compute_spie_simplified(I_theta);
[S0, S1, S2, DoLP_img, AoP_img] = compute_spie_simplified(I_theta);

if plotComparison; plot_pol_parameters_comparison(S0, S1, S2, DoLP_img, AoP_img, selectedType, 'SPIE Simple', outputFolder); end

%% Display
plot_polarization_spatially(DoLP_img, 'jet', ['SPIE Simplified Fourier Spatial Distribution of DoLP for ' selectedType], outputFolder);
plot_polarization_spatially(AoP_img, 'hsv', ['SPIE Simplified Fourier Spatial Distribution of AoLP for ' selectedType], outputFolder);

%% Save results
spie_simple_fourier_data_struct = struct('SpieSimple_Fourier_DoLP_map', DoLP_img, 'SpieSimple_Fourier_AoLP_map', AoP_img, 'Stokes_S0', S0, 'Stokes_S1', S1, 'Stokes_S2', S2, 'wavelenghts', wavelengths);

% Save the full struct to file
save_analysis_data(spie_simple_fourier_data_struct, 'SpieSimple_Fourier_Stokes_Results', selectedType, outputFolder);

%% Compare methods
compare_methods_DoLP(DoLP_wg, 'Standard', DoLP_fourier_img, 'Fourier', DoLP_img, 'SPIE simplified');

compare_methods_AoLP_hsi(AoLP_wg, 'Standard', AoLP_fourier_img, 'Fourier', AoP_img, 'SPIE simplified'); 

fprintf('\nStandard: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(DoLP_wg(:)), max(DoLP_wg(:)), (min(AoLP_wg(:))), (max(AoLP_wg(:))))
fprintf('\nFourier: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(DoLP_fourier_img(:)), max(DoLP_fourier_img(:)), (min(AoLP_fourier_img(:))), (max(AoLP_fourier_img(:))))
fprintf('\n SPIE Fourier: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(DoLP(:)), max(DoLP(:)), (min(AoP_deg(:))), (max(AoP_deg(:))))
fprintf('\n SPIE Fourier simplified: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(DoLP_img(:)), max(DoLP_img(:)), (min(AoP_img(:))), (max(AoP_img(:))))

end


