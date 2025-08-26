%% MAIN SCRIPT FOR HSI + POLARIZATION FUSION and Overlay Analysis - FINAL VERSION

clc; clear; close all;

%%
check_calibration = true;
offsetY = -75; % empirical correction

%% Get auxiliary matlab functions previously created for hsi and polarization data analysis, as well as auxiliary functions for this script
addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'pol'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi'));

%% Load dataset info
baseDir  = "D:\afili\Transferências\Tese - files\*";
basePath = fileparts(baseDir);
darkPath = 'dark.raw';

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

mainDatasetName = sprintf('%s_%s', selectedType, selectedDay);

%% Select methods
hsiMethod = questdlg('Select HSI method:', 'Method Selection (HSI)', 'Standard', 'Fourier', 'SPIE Simplified', 'Standard');

polMethod = questdlg('Select POL method:', 'Method Selection (POL)', 'Standard', '2nd Order Fourier', '4th Order Fourier', 'Standard');

fusionChoice = questdlg('Select fusion setup:', 'Fusion Options', 'HSI with WG + POL', 'HSI without WG + POL', 'HSI with WG + POL');

%% Prepare calibration images
processPolTestData(basePath, selectedType, selectedDay);
load('processedPolRawFiles.mat', 'calibPolData', 'woWGCalibData');

fprintf("Loading Checkerboard images...\n");
[hsi_rgb_u8, pol_rgb] = prepareCalibrationImages(basePath, calibHsData, calibPolData, selectedType, selectedDay, mainDatasetName);
%size(hsi_rgb_u8) % 921 x 1776 x 3
%size(pol_rgb) % 1024 x 1224 x 3

%% Run analysis based on previously obtained results
switch fusionChoice
    case 'HSI with WG + POL'
        [DoLP_HSI, AoLP_HSI, wavelengths] = loadHsiAnalysis(hsiMethod, selectedType, mainDatasetName);
        [DoLP_POL, AoLP_POL] = loadPolAnalysis(polMethod, mainDatasetName);
        %[DoLP_HSI_calib, ~, ~, resizeSize] = loadHsiAnalysis(hsiMethod, selectedType,  [selectedType selectedDay 'juncalib' ]);
        %[DoLP_POL_calib, ~] = loadPolAnalysis(polMethod,[selectedType '_calib_' selectedDay]);
        
        %% Calibration: Get homography from checkerboard
        %tform = calibrateHsiPolDataFromCheckerboard(DoLP_HSI_calib, DoLP_POL_calib, [selectedType '_calib_' selectedDay], size(DoLP_HSI), size(DoLP_POL), check_calibration);
        tform = calibrateHsiPolDataFromCheckerboard(hsi_rgb_u8, pol_rgb, [selectedType '_calib_' selectedDay], size(DoLP_HSI), size(DoLP_POL), check_calibration);

    case 'HSI without WG + POL'
        % Only reflectance fusion possible
        [reflectance_HSI, ~, wavelengths] = loadHsiAnalysis('Reflectance', selectedType, [selectedType selectedDay 'junwoWGpol']);
        [DoLP_POL, AoLP_POL] = loadPolAnalysis(polMethod, [selectedType '_woWG_' selectedDay]);
        %[reflectance_HSI_calib, ~,~, resizeSize] = loadHsiAnalysis('Reflectance', selectedType, [selectedType selectedDay 'junwoWGpolcalib' ]);
        %[DoLP_POL_calib, ~] = loadPolAnalysis(polMethod,[selectedType '_woWG_calib_' selectedDay]);
      
        %% Calibration: Get homography from checkerboard
        %tform = calibrateHsiPolDataFromCheckerboard(reflectance_HSI_calib, DoLP_POL_calib, [selectedType '_woWG_calib_' selectedDay], resizeSizeHsi, size(DoLP_POL), check_calibration);
        tform = calibrateHsiPolDataFromCheckerboard(hsi_rgb_u8, pol_rgb, [selectedType '_woWG_calib_' selectedDay], size(reflectance_HSI), size(DoLP_POL), check_calibration);
end

%% Calibration: Get homography from checkerboard
fprintf('Homography transformation matrix (tform.T) for "%s" dataset:\n', mainDatasetName);
disp(tform.T);

figure('Name','Comparison - checkerboard'); %before misalignment correction
subplot(1,3,1); imshow(mean(hsi_rgb_u8,3), []); colormap('jet'); colorbar; title('DoLP (HSI)');
subplot(1,3,2); imshow(mean(pol_rgb,3), []); colormap('jet'); colorbar; title('DoLP (POL)');
subplot(1,3,3); imshowpair(mean(imwarp(pol_rgb, tform, 'OutputView', imref2d(size(hsi_rgb_u8))),3),mean(hsi_rgb_u8,3), 'Diff'); title('Overlay checkerboard');
%% Alignment based on tform
if exist('DoLP_HSI','var')
    tform.T(3,2) = tform.T(3,2) + offsetY; % shift in Y

    aligned_DoLP_POL = imwarp(DoLP_POL, tform, 'OutputView', imref2d(size(DoLP_HSI)));
    aligned_AoLP_POL = imwarp(AoLP_POL, tform, 'OutputView', imref2d(size(AoLP_HSI))); % in rad
else
    aligned_DoLP_POL = imwarp(DoLP_POL, tform, 'OutputView', imref2d(size(reflectance_HSI)));
    aligned_AoLP_POL = imwarp(AoLP_POL, tform, 'OutputView', imref2d(size(reflectance_HSI))); % in rad
end

%% Visualization
if exist('DoLP_HSI','var')
    % Quick test: overlay DoLP with RGB reference to see if different processing results in disalignments
    hsi_rgb_u8 = imresize(hsi_rgb_u8, [size(DoLP_HSI,1), size(DoLP_HSI,2)]); %resize image based on reference
    
    figure('Name','Check HSI alignment');
    DoLP_norm = imbinarize(mean(DoLP_HSI,3), 'adaptive'); % Normalize DoLP for visibility
    DoLP_norm = imtranslate(DoLP_norm, [0, - offsetY]); % shift up by 80 px to account for vertical misalignment
    subplot(1,2,1); imshowpair(DoLP_norm, mean(hsi_rgb_u8,3),'falsecolor'); title('HSI RGB vs DoLP (HSI) - falsecolor overlay');
    subplot(1,2,2); imshowpair(mean(hsi_rgb_u8,3), DoLP_norm, 'diff'); title('HSI RGB vs DoLP (HSI) - diff overlay');
    
    figure('Name','Check Pol alignment');
    % Normalize DoLP for visibility
    DoLP_norm = imbinarize(mean(DoLP_POL,3), 'adaptive');
    subplot(1,2,1); imshowpair(DoLP_norm, mean(pol_rgb,3),'falsecolor'); title('Pol RGB vs DoLP (Pol) - falsecolor overlay');
    subplot(1,2,2); imshowpair(mean(pol_rgb,3), DoLP_norm, 'diff'); title('Pol RGB vs DoLP (Pol) - diff overlay');

    %% DoLP 
    figure('Name','Comparison - DoLP');
    subplot(1,3,1); imshow(mean(DoLP_HSI,3), []); colormap('jet'); colorbar; title('DoLP (HSI)'); axis image; caxis([0 1]);
    subplot(1,3,2); imshow(mean(aligned_DoLP_POL,3), []); colormap('jet'); colorbar; title('DoLP (POL)'); axis image; caxis([0 1]);
    subplot(1,3,3); imshowpair(mean(DoLP_HSI,3), mean(aligned_DoLP_POL,3), 'Diff'); title('Overlay DoLP'); axis image; %caxis([0 1]);

    figure('Name','HSI POL Fusion - DoLP');
    % Use false color overlay with transparency
    imshow(DoLP_HSI, [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(aligned_DoLP_POL, [0 1]);
    set(h, 'AlphaData', 0.65);  % Semi-transparent overlay
    title('Overlay DoLP'); axis image; %caxis([0 1]);
    
    %% AoLP
    figure; imshowpair(AoLP_HSI, rad2deg(mean(aligned_AoLP_POL,3)), 'falsecolor');
    title('AoLP (HSI) vs AoLP (POL)');

    figure('Name','Comparison - AoLP');
    subplot(1,3,1); imshow(mean(AoLP_HSI,3), []); colormap('hsv'); colorbar; title('AoLP (HSI)'); axis image; caxis([-90 90]);
    subplot(1,3,2); imshow(rad2deg(mean(aligned_AoLP_POL,3)), []); colormap('hsv'); colorbar; title('AoLP (POL)'); axis image; %caxis([-90 90]);
    subplot(1,3,3); imshowpair(mean(AoLP_HSI,3), rad2deg(mean(aligned_AoLP_POL,3)), 'Diff'); title('Overlay AoLP'); axis image; %caxis([-90 90]);

    %% Smoothed and brightened images

    bw_DoLP = uint8(255 * mat2gray(aligned_DoLP_POL, [prctile(aligned_DoLP_POL(:), 1), prctile(aligned_DoLP_POL(:), 99)]));
    bw_AoLP = uint8(255 * mat2gray(rad2deg(aligned_AoLP_POL), [prctile(rad2deg(aligned_DoLP_POL(:)), 1), prctile(rad2deg(aligned_DoLP_POL(:)), 99)]));  % [0, 255]
     
    bw_DoLP_smooth = imfilter(bw_DoLP, fspecial('average', [2 2])); % the larger the kernel, the blurrier the image
    bw_AoLP_smooth = imfilter(bw_AoLP, fspecial('average', [3 3]));
    
    DoLP_HSI_smooth = imfilter(DoLP_HSI, fspecial('average', [2 2]));
    AoLP_HSI_smooth = imfilter(AoLP_HSI, fspecial('average', [2 2]));
    
    % Show results
    fig = figure('Name', ['Smoothed DoLP & AoLP (Grayscale) for ' mainDatasetName]);
    
    subplot(1,2,1);
    imshow(mean(DoLP_HSI_smooth,3), [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(bw_DoLP_smooth, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('DoLP (BW + Brightened + Smoothed)');
    
    subplot(1,2,2);
    imshow(mean(AoLP_HSI_smooth,3), []); colormap(jet); colorbar;axis image; hold on;
    h = imshow(bw_AoLP_smooth, [-90 90]); 
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('AoLP (BW + Brightened + Smoothed)');

else
    disp('Only POL vs Reflectance fusion available (HSI wo/ WG)');
    
    %% DoLP
    figure; imshowpair(reflectance_HSI, mean(aligned_DoLP_POL,3), 'falsecolor');
    title('Reflectance (HSI wo/ WG) vs DoLP (POL)');

    reflectance_adj = adapthisteq(mean(reflectance_HSI,3)); %add contrast
    
    figure('Name','Comparison - Reflectance and DoLP');
    subplot(1,3,1); imshow(reflectance_adj, []); colormap('jet'); colorbar; title('Reflectance (HSI)'); axis image; caxis([0 1]);
    subplot(1,3,2); imshow(mean(aligned_DoLP_POL,3), []); colormap('jet'); colorbar; title('DoLP (POL)'); axis image; caxis([0 1]);
    subplot(1,3,3); imshowpair(reflectance_adj, mean(aligned_DoLP_POL,3), 'Diff'); title('Reflectance with DoLP Overlay'); axis image; 

    %% AoLP
    figure; imshowpair(reflectance_HSI, mean(aligned_AoLP_POL,3), 'falsecolor');
    title('Reflectance (HSI wo/ WG) vs AoLP (POL)');

    figure('Name','Comparison - Reflectance and AoLP');
    subplot(1,3,1); imshow(reflectance_adj, []); colormap('hsv'); colorbar; title('Reflectance (HSI)'); axis image; %caxis([-90 90]);
    subplot(1,3,2); imshow(mean(aligned_AoLP_POL,3), []); colormap('hsv'); colorbar; title('AoLP (POL)'); axis image; %caxis([-90 90]);
    subplot(1,3,3); imshowpair(reflectance_adj, mean(aligned_AoLP_POL,3), 'Diff'); title('Reflectance with AoLP Overlay'); axis image; %caxis([-90 90]);


    %% Smoothed and brightened images
    
    bw_DoLP = uint8(255 * mat2gray(aligned_DoLP_POL, [prctile(aligned_DoLP_POL(:), 1), prctile(aligned_DoLP_POL(:), 99)]));
    bw_AoLP = uint8(255 * mat2gray(rad2deg(aligned_AoLP_POL), [prctile(rad2deg(aligned_DoLP_POL(:)), 1), prctile(rad2deg(aligned_DoLP_POL(:)), 99)]));  % [0, 255]
     
    bw_DoLP_smooth = imfilter(bw_DoLP, fspecial('average', [2 2])); % the larger the kernel, the blurrier the image
    bw_AoLP_smooth = imfilter(bw_AoLP, fspecial('average', [3 3]));
    
    % Show results
    fig = figure('Name', ['Smoothed DoLP & AoLP (Grayscale) for ' mainDatasetName]);
    
    subplot(1,2,1);
    imshow(reflectance_adj, [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(bw_DoLP_smooth, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('DoLP (BW + Brightened + Smoothed)');
    
    subplot(1,2,2);
    imshow(reflectance_adj, [0 1]); colormap(jet); colorbar;axis image; hold on;
    h = imshow(bw_AoLP_smooth, [-90 90]); 
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('AoLP (BW + Brightened + Smoothed)');
end



