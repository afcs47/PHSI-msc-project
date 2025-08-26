% HSI + Polarization Fusion and Overlay Analysis
clc; clear; close all;

%% Get auxiliary matlab functions previously created for HSI and polarization data analysis
addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'pol'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi'));

%% Define paths and dataset
baseDir = "D:\afili\Transferências\Tese - files\*";
basePath = fileparts(baseDir);
dark_path = 'dark.raw';

processHsTestData(baseDir);
load('allHsTestData.mat', 'stokesData', 'woWGpolData', 'fourierData', 'calibHsData');

%% Choose sample and day
sampleTypes = unique(stokesData.SampleName);
[selectedTypeIndex, tf] = listdlg('PromptString','Select sample type:', 'SelectionMode','single', 'ListString', sampleTypes);
if ~tf, return; end
selectedType = sampleTypes{selectedTypeIndex};

% Filter stokesData by selected type by comparing strings/names
filteredData = stokesData(strcmp(stokesData.SampleName, selectedType), :);

dayList = unique(filteredData.Day);
if length(dayList) > 1
    [selectedDayIndex, tf] = listdlg('PromptString','Select day:', 'SelectionMode','single', 'ListString', dayList);
    if ~tf, return; end
    selectedDay = dayList{selectedDayIndex};
else
    selectedDay = dayList{1};
end

%% Run HSI analysis to compute DoLP and AoLP
% Try loading precomputed HSI data
%hsiDataFile = dir(fullfile(['hsi results+figures\' selectedType '_' selectedDay], sprintf('SpieSimple_Fourier_Stokes_Results_%s_*.mat', selectedType)));
hsiDataFile = dir(fullfile(['hsi results+figures\' selectedType '_' selectedDay], sprintf('Standard_Stokes_Results_%s_*.mat', selectedType))); %Change file name according to the chosen computation method
%if ~isempty(hsiDataFile)
    %load(fullfile('hsi results+figures', hsiDataFile(end).name), 'SpieSimple_Fourier_DoLP_map', 'SpieSimple_Fourier_AoLP_map', 'Stokes_S0', 'Stokes_S1', 'Stokes_S2', 'wavelengths');
    %DoLP_HSI = SpieSimple_Fourier_DoLP_map;
    %AoLP_HSI = SpieSimple_Fourier_AoLP_map;

if isempty(hsiDataFile); error('Error: HSI data file for %s - %s not found! Please run the respective script first or check directories.', datasetname); end %comment when uncommenting other lines after if ~isempty
    load(fullfile(['hsi results+figures\' selectedType '_' selectedDay], hsiDataFile(end).name), 'Standard_DoLP_map', 'Standard_AoLP_map', 'Stokes_S0', 'Stokes_S1', 'Stokes_S2', 'wavelengths');
    DoLP_HSI = Standard_DoLP_map;
    AoLP_HSI = Standard_AoLP_map;

% else % Compute if not available
%     % Filter Fourier data by selected sample and day
%     filteredFourier = fourierData(strcmp(fourierData.SampleName, selectedType) & strcmp(fourierData.Day, selectedDay), :);
%     anglesDeg = cellfun(@str2double, filteredFourier.WGpolAngle);
%     [sortedAngles, sortIdx] = sort(anglesDeg);
%     filteredFourier = filteredFourier(sortIdx, :);
%     
%     % Compute I_theta cube
%     mean_ref_struct = struct();
%     for i = 1:numel(sortedAngles)
%         fpath = fullfile(basePath, filteredFourier.FolderName{i}, 'capture\');
%         fname = filteredFourier.FolderName{i};
%         [Data, White, Dark, ~] = read_data(fpath, fname);
%         HS_calibrated = apply_calibration(Data, White, Dark);
%         hsfiltered = apply_sg_filter(HS_calibrated);
%         mean_ref = mean(hsfiltered, 3);
%         fieldName = sprintf('angle_%d', sortedAngles(i));
%         mean_ref_struct.(fieldName) = mean_ref;
%     end
%     mean_ref_struct = resize_reflectances(mean_ref_struct);
%     fields = fieldnames(mean_ref_struct);
%     numAngles = numel(fields);
%     [rows, cols] = size(mean_ref_struct.(fields{1}));
%     I_theta = zeros(rows, cols, numAngles);
%     for i = 1:numAngles
%         I_theta(:, :, i) = mean_ref_struct.(fields{i});
%     end
%     
%     % Compute simplified SPIE method DoLP and AoLP
%     [DoLP_HSI, AoLP_HSI] = compute_spie_simplified(I_theta);
% end

%% Run Pol analysis to compute DoLP and AoLP
% Try loading polarization data
polDataFile = fullfile(['pol results+figures\' selectedType '_' selectedDay 'jun.raw'], sprintf('%s_%sjun.raw_polarization_data.mat', selectedType, selectedDay)); %Change file name according to chosen computation method
%if exist(polDataFile, 'file')
%     load(polDataFile, 'DoLP_2nd', 'AoLP_2nd');
%     DoLP_POL = DoLP_2nd;
%     AoLP_POL = AoLP_2nd;

if ~isfile(polDataFile); error('Error: Polarization data file for %s - %s not found! Please run the respective script first or check directories.', datasetname); end  %comment when uncommenting other lines after if ~isempty
    load(polDataFile, 'DoLP', 'AoLP');
    DoLP_POL = DoLP;
    AoLP_POL = AoLP;

% else
%     % Select respective polarization data
      processPolTestData(basePath, selectedType, selectedDay);
      load('processedPolRawFiles.mat', 'standardData', 'calibPolData', 'woWGData', 'woWGCalibData')
%     
%     % When there's more than one raw file per sample, let user select one
%     versionList = standardData.SampleVersion;
%     if length(versionList)>1
%         [selectedVersionIndex, tf] = listdlg('PromptString','Select sample version:', 'SelectionMode','single', 'ListString', versionList);
%         if ~tf, return; end
%         standardData = standardData(strcmp(standardData.SampleVersion, versionList{selectedVersionIndex}), :);
%     end
%     
%     input_pol_raw = fullfile(basePath, "pol data/", standardData.FileName);
%     
%     % Run polarization analysis
%     Z = load_image(input_pol_raw, '', 0, 0, dark_path);
%     Z = imrotate(Z, 180);
%     [proc_90, proc_45, proc_135, proc_0] = demosaic_polarization_image(Z, input_pol_raw, '', 0, 0);
%     [DoLP_POL, AoLP_POL] = calculate_polarization(proc_90, proc_45, proc_135, proc_0);
%     DoLP_POL_mean = mean(DoLP_POL, 3);
%     AoLP_POL_mean = mean(AoLP_POL, 3);
% end
%%
figure;
imshowpair(mean(DoLP_POL,3), mean(DoLP_HSI,3),'Diff') %shows initial disalignment

%% Calibration: Get homography from checkerboard
selectedCalibData = calibHsData(strcmp(calibHsData.SampleName, selectedType) & strcmp(calibHsData.Day, selectedDay), :);
fname = selectedCalibData.FolderName{1};
fpath = fullfile(basePath, fname, 'capture\');
datasetname = sprintf('%s_%s', selectedType, selectedDay);
calib_pol_raw = fullfile(basePath, "pol data/", calibPolData.FileName);
%%
tform = calibrateHsiPolData(fpath, fname, calib_pol_raw, datasetname, size(DoLP_HSI), size(DoLP_POL), true);

fprintf('Homography transformation matrix (tform.T) for "%s" dataset:\n', datasetname);
disp(tform.T);


%% Align POL images to HSI frame
aligned_DoLP_POL = imwarp(DoLP_POL, tform, 'OutputView', imref2d(size(DoLP_HSI)));
aligned_AoLP_POL = rad2deg(imwarp(AoLP_POL, tform, 'OutputView', imref2d(size(AoLP_HSI))));
rgb_aligned_AoLP_POL = map_aolp_to_rgb(mean(aligned_AoLP_POL,3), -90, 90);

% Overlay and display of DoLP
figure('Name','Comparison - DoLP');
subplot(1,3,1); imshow(mean(DoLP_HSI,3), []); colormap('jet'); colorbar; title('DoLP (HSI)'); axis image; caxis([0 1]);
subplot(1,3,2); imshow(mean(aligned_DoLP_POL,3)); colormap('jet'); colorbar; title('DoLP (POL)'); axis image; caxis([0 1]);
subplot(1,3,3); imshowpair(mean(DoLP_HSI,3), mean(aligned_DoLP_POL,3), 'Diff'); colormap('jet'); colorbar; title('Overlay DoLP'); axis image; %caxis([0 1]);
% Use false color overlay with transparency
imshow(DoLP_HSI, [0 1]); colormap(turbo); colorbar; axis image; hold on;
h = imshow(aligned_DoLP_POL, [0 1]);
set(h, 'AlphaData', 0.5);  % Semi-transparent overlay

%% AoLP
figure('Name','Comparison - AoLP');
subplot(1,3,1); imshow(mean(AoLP_HSI,3), []); colormap('hsv'); colorbar; title('AoLP (HSI)'); axis image; %caxis([-90 90]);
subplot(1,3,2); imshow(rgb_aligned_AoLP_POL); colormap('hsv'); colorbar; title('AoLP (POL)'); axis image; %caxis([-90 90]);
subplot(1,3,3); imshowpair(mean(AoLP_HSI,3), aligned_AoLP_POL, 'Diff'); colormap('hsv'); colorbar; title('Overlay AoLP'); axis image; %caxis([-90 90]);

%% Save fused data
save(['Fused_DoLP_AoLP_' datasetname '.mat'], 'DoLP_HSI', 'AoLP_HSI', 'aligned_DoLP_POL', 'aligned_AoLP_POL');
