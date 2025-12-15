%% MAIN SCRIPT FOR HSI + POLARIZATION FUSION and Overlay Analysis - FINAL VERSION

clc; clear; close all;

%%
check_calibration = true;
%offsetY = -75; % empirical correction

%% Get auxiliary matlab functions previously created for hsi and polarization data analysis, as well as auxiliary functions for this script
addpath(fullfile(fileparts(mfilename('fullpath')), 'hsi'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'pol'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'phsi'));

%% Start logging all console output
% Create output folder
outputFolderOG = uigetdir(pwd, 'Select Output Folder'); % Ask where to save results
if outputFolderOG == 0
    disp('No output folder selected. Exiting.');
    return
end

%% Load dataset info
baseDir  = "D:\afili\Transferências\Tese - files\*";
basePath = fileparts(baseDir);
darkPath = 'dark.raw';

process_hs_testData(baseDir);
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

%% Create text file to save command window output
mainDatasetName = sprintf('%s_%s', selectedType, selectedDay); %woWG data has a different day from hsi 

outputFolder = fullfile(outputFolderOG, selectedType);
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end % Create the subfolder based on the filename variable if it doesn't exist

logFile = fullfile(outputFolder, ['runlog_' + string(datetime("now","Format","yyyyMMdd_HHmmss")) + '.txt']);
diary(logFile);
diary on;
disp('=== Command Window Log Started ===');

disp(['Dataset: ' selectedType]);

%% Select methods
fusionChoice = questdlg('Select fusion setup:', 'Fusion Options', 'HSI with WG + POL', 'HSI without WG + POL', 'HSI with WG + POL');
disp(['Fusion choice: ' fusionChoice]);

polMethod = questdlg('Select POL method:', 'Method Selection (POL)', 'Standard', '2nd Order Fourier', '4th Order Fourier', 'Standard');
disp(['POL method: ' polMethod]);

% HSI method selection: allow many processing choices for rotating WG
hsiMethod = '';
if strcmp(fusionChoice, 'HSI with WG + POL')
    hsiMethod = questdlg('Select HSI method:', 'Method Selection (HSI)', 'Standard', 'FFT', 'LSQ', 'Standard');
    % last item returned as default if user closes, ensure variable
    if isempty(hsiMethod)
        hsiMethod = 'Standard';
    end
    disp(['HSI method: ' hsiMethod]);
end

%% Prepare calibration images
process_pol_testData(basePath, selectedType, selectedDay);
load('processedPolRawFiles.mat', 'calibPolData', 'woWGCalibData');

fprintf("Loading Checkerboard images...\n");
[hsi_rgb_u8, pol_rgb] = prepare_calibration_images(basePath, calibHsData, calibPolData, selectedType, selectedDay, mainDatasetName);
%size(hsi_rgb_u8) % 921 x 1776 x 3
%size(pol_rgb) % 1024 x 1224 x 3

%% Run analysis based on previously obtained results
switch fusionChoice
    case 'HSI with WG + POL'
        %hsiMethod = questdlg('Select HSI method:', 'Method Selection (HSI)', 'Standard', 'FFT', 'LSQ', 'Standard');
        %disp(['HSI method: ' hsiMethod]);

        [DoLP_HSI, AoLP_HSI] = load_hsi_analysis(hsiMethod, selectedType, mainDatasetName);
        [DoLP_POL, AoLP_POL] = load_pol_analysis(polMethod, mainDatasetName);
        [cube, wavelengths] = load_hsi_analysis('Reflectance', selectedType, mainDatasetName, 'Spectral_reflectances');

        %% Calibration: Get homography from checkerboard
        tform = calibrate_hsi_polData_from_checkerboard(hsi_rgb_u8, pol_rgb, [selectedType '_calib_' selectedDay], size(DoLP_HSI), size(DoLP_POL), check_calibration);
        
    case 'HSI without WG + POL'
        % Only reflectance fusion possible
        [reflectance_HSI, ~] = load_hsi_analysis('Reflectance', selectedType, [selectedType '_' woWGCalibData.Day{1} '_woWG'], 'Spatial_reflectances'); %selectedDay of sampled data = day of calibration data
        [DoLP_POL, AoLP_POL] = load_pol_analysis(polMethod, [selectedType '_woWG_' woWGCalibData.Day{1}]);
        [cube, wavelengths] = load_hsi_analysis('Reflectance', selectedType, [selectedType '_' woWGCalibData.Day{1} '_woWG'], 'Spectral_reflectances');
        reflectance_HSI = reflectance_HSI.sample_01;
        cube = cube.sample_01;
        %% Calibration: Get homography from checkerboard
        tform = calibrate_hsi_polData_from_checkerboard(hsi_rgb_u8, pol_rgb, [selectedType '_woWG_calib_' selectedDay], size(reflectance_HSI), size(DoLP_POL), check_calibration);
end

% offsetY = estimate_offsetY(mean(hsi_rgb_u8,3), mean(pol_rgb,3));
% fprintf('User-estimated offsetY = %d pixels\n', offsetY);


%% Calibration: Get homography from checkerboard
fprintf('Homography transformation matrix (tform.T) for "%s" dataset:\n', mainDatasetName);
disp(tform.T);

figure('Name','Comparison - checkerboard'); %before misalignment correction
subplot(1,3,1); imshow(mean(hsi_rgb_u8,3), []); title('Checkerboard (HSI)');
subplot(1,3,2); imshow(mean(pol_rgb,3), []); title('Checkerboard (POL)');
subplot(1,3,3); imshowpair(mean(imwarp(pol_rgb, tform, 'OutputView', imref2d(size(hsi_rgb_u8))),3),mean(hsi_rgb_u8,3), 'Diff'); title('Overlay checkerboard');


%% Alignment based on tform (without accounting for y-offset)
if exist('DoLP_HSI','var')
    %tform.T(3,2) = tform.T(3,2) + offsetY; % shift in Y

    aligned_DoLP_POL = imwarp(DoLP_POL, tform, 'OutputView', imref2d(size(DoLP_HSI)));
    aligned_AoLP_POL = imwarp(AoLP_POL, tform, 'OutputView', imref2d(size(AoLP_HSI))); % in rad
else
    aligned_DoLP_POL = imwarp(DoLP_POL, tform, 'OutputView', imref2d(size(reflectance_HSI)));
    aligned_AoLP_POL = imwarp(AoLP_POL, tform, 'OutputView', imref2d(size(reflectance_HSI))); % in rad
end


%% Interactive offset adjustment - Correction of the DoLP/AoLP fusion vertical offset, not accounted by homography matrix (replaces fixed offsetY)
if exist('DoLP_HSI','var')
    disp('Estimating offsetY between DoLP (HSI) and warped DoLP (POL)...');
    offsetY = estimate_offsetY(DoLP_HSI, aligned_DoLP_POL);
    fprintf('User-estimated offsetY = %d pixels\n', offsetY);

    % Apply the additional Y-shift to the already warped POL data
    aligned_DoLP_POL = imtranslate(aligned_DoLP_POL, [0 offsetY]);
    aligned_AoLP_POL = imtranslate(aligned_AoLP_POL, [0 offsetY]);
else
    disp('Estimating offsetY between Reflectance (HSI) and warped DoLP (POL)...');
    offsetY = estimate_offsetY(reflectance_HSI, aligned_DoLP_POL);
    fprintf('User-estimated offsetY = %d pixels\n', offsetY);

    aligned_DoLP_POL = imtranslate(aligned_DoLP_POL, [0 offsetY]);
    aligned_AoLP_POL = imtranslate(aligned_AoLP_POL, [0 offsetY]);
end



%% Visualization
if exist('DoLP_HSI','var')
    % Quick test: overlay DoLP with RGB reference to see if different processing results in disalignments
    hsi_rgb_u8 = imresize(hsi_rgb_u8, [size(DoLP_HSI,1), size(DoLP_HSI,2)]); %resize image based on reference
    
    figure('Name','Check HSI alignment');
    DoLP_norm = imbinarize(mean(DoLP_HSI,3), 'adaptive'); % Normalize DoLP for visibility
    DoLP_norm = imtranslate(DoLP_norm, [0, - offsetY]); % shift up to account for vertical misalignment
    subplot(1,2,1); imshowpair(DoLP_norm, mean(hsi_rgb_u8,3),'falsecolor'); title('HSI RGB vs DoLP (HSI) - falsecolor overlay');
    subplot(1,2,2); imshowpair(mean(hsi_rgb_u8,3), DoLP_norm, 'diff'); title('HSI RGB vs DoLP (HSI) - diff overlay');
    
    figure('Name','Check Pol alignment');
    % Normalize DoLP for visibility
    DoLP_norm = imbinarize(mean(DoLP_POL,3), 'adaptive');
    subplot(1,2,1); imshowpair(DoLP_norm, mean(pol_rgb,3),'falsecolor'); title('Pol RGB vs DoLP (Pol) - falsecolor overlay');
    subplot(1,2,2); imshowpair(mean(pol_rgb,3), DoLP_norm, 'diff'); title('Pol RGB vs DoLP (Pol) - diff overlay');

    %% DoLP 
    figure('Name','Comparison - DoLP');
    subplot(1,2,1); imshow(mean(DoLP_HSI,3), []); colormap('jet'); colorbar; title('DoLP (HSI)'); axis image; caxis([0 1]);
    subplot(1,2,2); imshow(mean(aligned_DoLP_POL,3), []); colormap('jet'); colorbar; title('DoLP (POL)'); axis image; caxis([0 1]);
    title("HSI vs POL DoLP");

    
    figure('Name','HSI POL Fusion - DoLP');
    % Use false color overlay with transparency
    imshow(DoLP_HSI, [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(aligned_DoLP_POL, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('Overlay DoLP'); axis image; %caxis([0 1]);


    %% Brightened images
    bw_DoLP = uint8(255 * mat2gray(aligned_DoLP_POL, [prctile(aligned_DoLP_POL(:), 1), prctile(aligned_DoLP_POL(:), 99)]));
    
    % Show results
    figure('Name', ['Brightened DoLP for ' mainDatasetName]);
    imshow(mean(DoLP_HSI,3), [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(bw_DoLP, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('DoLP (Brightened)');

    %%
    DoLP_hsi_rgb = repmat(DoLP_HSI, 1, 1, 3);  % convert grayscale to RGB
    bw_DoLP = bw_DoLP(:,:,1);   % ensure single plane
    dolp_rgb = repmat(bw_DoLP(:,:,1), 1, 1, 3); % grayscale to RGB

    figure('Name', ['Brightened DoLP (Grayscale) for ' mainDatasetName]);
    imshow(mean(DoLP_hsi_rgb,3), [0 1]); colormap("gray"); colorbar; axis image; hold on;
    h = imshow(dolp_rgb, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('DoLP (BW + Brightened)');


    %% Convert DoLP maps to GRAYSCALE (not jet)
    clims = [0 1];
    
    % Reduce to 2D DoLP maps
    DoLP_HSI_gray = mat2gray(mean(DoLP_HSI, 3), clims);
    DoLP_POL_gray = mat2gray(mean(aligned_DoLP_POL, 3), clims);
    
    % Convert to grayscale RGB (MxNx3)
    DoLP_HSI_rgb = repmat(DoLP_HSI_gray, 1, 1, 3);
    DoLP_POL_rgb = repmat(DoLP_POL_gray, 1, 1, 3);
    
    % Side-by-side comparison 
    figure('Name','Comparison - DoLP (Grayscale)');
    montage({DoLP_HSI_rgb, DoLP_POL_rgb}, 'Size', [1 2]);
    title("HSI vs POL DoLP (Grayscale)");
    
    % Grayscale fusion 
    figure('Name','Fused DoLP (Grayscale)');
    imshow(DoLP_HSI_rgb); hold on;
    h = imshow(DoLP_POL_rgb);
    set(h, 'AlphaData', 0.5);   % transparency
    title('DoLP HSI + DoLP POL fused (Grayscale)');
    axis image;

    %% Convert DoLP maps to RGB using same colormap + limits
    cmap = jet(256);% choose the same colormap
    clims = [0 1]; % consistent color range for both modalities
    DoLP_HSI = mean(DoLP_HSI, 3);
    aligned_DoLP_POL = mean(aligned_DoLP_POL, 3);
    
    DoLP_HSI_rgb = ind2rgb(gray2ind(mat2gray(DoLP_HSI, clims), 256), cmap);
    DoLP_POL_rgb = ind2rgb(gray2ind(mat2gray(aligned_DoLP_POL, clims), 256), cmap);
    
    figure('Name','Comparison - DoLP 2');
    montage({DoLP_HSI_rgb, DoLP_POL_rgb}, 'Size', [1 2]); %side-by-side without colorbar
    title("HSI vs POL DoLP (same color scale)");

    fusedFig = figure('Name','Fused DoLP');
    imshow(DoLP_HSI_rgb); hold on;
    h = imshow(DoLP_POL_rgb);
    transpValue = 0.5;
    set(h, 'AlphaData', transpValue); % adjust transparency
    title(sprintf('DoLP HSI + POL fused (consistent colormap, transparency @ %g)',transpValue));
    axis image;
    hold on;
    
%     %% AoLP
%     figure('Name','HSI POL Fusion - AoLP'); imshowpair(AoLP_HSI, rad2deg(mean(aligned_AoLP_POL,3)), 'falsecolor'); % doesn't look great
%     title('AoLP (HSI) vs AoLP (POL)');
% 
%     figure('Name','Comparison - AoLP');
%     subplot(1,3,1); imshow(mean(AoLP_HSI,3), []); colormap('hsv'); colorbar; title('AoLP (HSI)'); axis image; caxis([-90 90]);
%     subplot(1,3,2); imshow(rad2deg(mean(aligned_AoLP_POL,3)), []); colormap('hsv'); colorbar; title('AoLP (POL)'); axis image; %caxis([-90 90]);
%     subplot(1,3,3); imshowpair(mean(AoLP_HSI,3), rad2deg(mean(aligned_AoLP_POL,3)), 'Diff'); title('Overlay AoLP'); axis image; %caxis([-90 90]);
% 
%     %% Convert AoLP maps to RGB - still doesn't result in a good fused image - harder since both aolp data is plotted from different functions
%     AoLP_HSI = mean(AoLP_HSI, 3); 
%     aligned_AoLP_POL = mean(aligned_AoLP_POL, 3);
%     
%     % Build HSV channels: hue = angle, saturation = 1, value = 1
%     AoLP_HSI_rgb = hsv2rgb(cat(3, AoLP_HSI, ones(size(AoLP_HSI)), ones(size(AoLP_HSI))));
%     AoLP_POL_rgb = hsv2rgb(cat(3, aligned_AoLP_POL, ones(size(aligned_AoLP_POL)), ones(size(aligned_AoLP_POL))));
%     
%     figure('Name','AoLP Fusion');
%     imshow(AoLP_HSI_rgb); hold on;
%     h = imshow(AoLP_POL_rgb);
%     set(h,'AlphaData', 0.35);
%     title('AoLP Fusion');
%     axis image;
  
else
    disp('Only POL vs Reflectance fusion available (HSI wo/ WG)');
    
    % Quick test: overlay DoLP with RGB reference to see if different processing results in disalignments
    hsi_rgb_u8 = imresize(hsi_rgb_u8, [size(reflectance_HSI,1), size(reflectance_HSI,2)]); %resize image based on reference
    
    figure('Name','Check HSI alignment');
    Reflec_norm = imbinarize(reflectance_HSI, 'adaptive');
    Reflec_norm = imtranslate(Reflec_norm, [0, - offsetY]);
    subplot(1,2,1); imshowpair(Reflec_norm, mean(hsi_rgb_u8,3),'falsecolor'); title('HSI RGB vs DoLP (HSI) - falsecolor overlay');
    subplot(1,2,2); imshowpair(mean(hsi_rgb_u8,3), Reflec_norm, 'diff'); title('HSI RGB vs DoLP (HSI) - diff overlay');
    
    figure('Name','Check Pol alignment');
    % Normalize DoLP for visibility
    DoLP_norm = imbinarize(mean(DoLP_POL,3), 'adaptive');
    subplot(1,2,1); imshowpair(DoLP_norm, mean(pol_rgb,3),'falsecolor'); title('Pol RGB vs DoLP (Pol) - falsecolor overlay');
    subplot(1,2,2); imshowpair(mean(pol_rgb,3), DoLP_norm, 'diff'); title('Pol RGB vs DoLP (Pol) - diff overlay');


    %% DoLP
    reflectance_adj = adapthisteq(reflectance_HSI); %adjust contrast to reflectance plot
    
    figure('Name','Comparison - Reflectance and DoLP');
    %subplot(1,3,1); imshow(mean(reflectance_HSI,3), []); colormap('jet'); colorbar; title('Reflectance (HSI)'); axis image; caxis([0 1]); %low contrast 
    subplot(1,2,1); imshow(reflectance_adj, []); colormap('jet'); colorbar; title('Reflectance (HSI)'); axis image; caxis([0 1]);
    subplot(1,2,2); imshow(mean(aligned_DoLP_POL,3), []); colormap('jet'); colorbar; title('DoLP (POL)'); axis image; caxis([0 1]);
    %subplot(1,3,3); imshowpair(reflectance_adj, mean(aligned_DoLP_POL,3), 'Diff'); title('Reflectance + DoLP'); axis image; 
    title("HSI Reflectance vs POL DoLP");

    figure('Name','HSI POL Fusion - Reflectance and DoLP');
    % Use false color overlay with transparency
    imshow(reflectance_adj, [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(aligned_DoLP_POL, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('Overlay Reflectance and DoLP'); axis image; %caxis([0 1]);

    %% Brightened images
    bw_DoLP = uint8(255 * mat2gray(aligned_DoLP_POL, [prctile(aligned_DoLP_POL(:), 1), prctile(aligned_DoLP_POL(:), 99)]));
    %bw_DoLP_smooth = imfilter(bw_DoLP, fspecial('average', [2 2])); % the larger the kernel, the blurrier the image
    % Show results
    figure('Name', ['Reflectance + DoLP for ' mainDatasetName]);
    imshow(mean(reflectance_adj,3), [0 1]); colormap(jet); colorbar; axis image; hold on;
    h = imshow(bw_DoLP, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('Reflectance and DoLP (Brightened)');

    %%
    refl_rgb = repmat(reflectance_adj, 1, 1, 3);  % convert grayscale to RGB
    bw_DoLP = bw_DoLP(:,:,1);   % ensure single plane
    dolp_rgb = repmat(bw_DoLP(:,:,1), 1, 1, 3); % grayscale to RGB

    figure('Name', ['Reflectance + DoLP (Grayscale) for ' mainDatasetName]);
    imshow(mean(refl_rgb,3), [0 1]); colormap("gray"); colorbar; axis image; hold on;
    h = imshow(dolp_rgb, [0 1]);
    set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
    title('Reflectance and DoLP (BW + Brightened)');

    %% Convert DoLP maps to GRAYSCALE (not jet)
    clims = [0 1];
    
    % Reduce to 2D DoLP maps
    Reflec_HSI_gray = mat2gray(reflectance_HSI, clims);
    DoLP_POL_gray = mat2gray(mean(aligned_DoLP_POL, 3), clims);
    
    % Convert to grayscale RGB (MxNx3)
    Reflec_HSI_rgb = repmat(Reflec_HSI_gray, 1, 1, 3);
    DoLP_POL_rgb = repmat(DoLP_POL_gray, 1, 1, 3);
    
    % Side-by-side comparison
    figure('Name','Comparison - DoLP (Grayscale)');
    montage({Reflec_HSI_rgb, DoLP_POL_rgb}, 'Size', [1 2]);
    title("HSI vs POL DoLP (Grayscale)");
    
    % Grayscale fusion 
    figure('Name','Fused DoLP (Grayscale)');
    imshow(Reflec_HSI_rgb); hold on;
    h = imshow(DoLP_POL_rgb);
    set(h, 'AlphaData', 0.5);   % transparency
    title('DoLP HSI + DoLP POL fused (Grayscale)');
    axis image;

    %% Convert DoLP maps to RGB using same colormap + limits
    cmap = jet(256);% choose the same colormap
    clims = [0 1]; % consistent color range for both modalities
    reflectance_adj = mean(reflectance_adj, 3);
    aligned_DoLP_POL = mean(aligned_DoLP_POL, 3);

    reflect_HSI_rgb = ind2rgb(gray2ind(mat2gray(reflectance_adj, clims), 256), cmap);
    DoLP_POL_rgb = ind2rgb(gray2ind(mat2gray(aligned_DoLP_POL, clims), 256), cmap);
    
    figure('Name','Comparison - Reflectance and DoLP 2');
    montage({reflect_HSI_rgb, DoLP_POL_rgb}, 'Size', [1 2]); %side-by-side without colorbar
    title("HSI Reflectance vs POL DoLP (same color scale)");

    fusedFig = figure('Name','Fused DoLP');
    imshow(reflect_HSI_rgb, [0 1]); hold on;
    h = imshow(DoLP_POL_rgb, [0 1]);
    transpValue = 0.35;
    set(h, 'AlphaData', transpValue); % adjust transparency
    title(sprintf('HSI Reflectance + DoLP POL fused (consistent colormap, transparency @ %g)',transpValue));
    axis image;
    hold on;

    
    %% AoLP
%     figure('Name','Comparison - Reflectance and AoLP');
%     subplot(1,3,1); imshow(reflectance_adj, []); colormap('hsv'); colorbar; title('Reflectance (HSI)'); axis image; %caxis([-90 90]);
%     subplot(1,3,2); imshow(mean(aligned_AoLP_POL,3), []); colormap('hsv'); colorbar; title('AoLP (POL)'); axis image; %caxis([-90 90]);
%     subplot(1,3,3); imshowpair(reflectance_adj, rad2deg(mean(aligned_AoLP_POL,3)), 'Diff'); title('Reflectance + AoLP'); axis image; %caxis([-90 90]);


%     bw_AoLP = uint8(255 * mat2gray(rad2deg(aligned_AoLP_POL), [prctile(rad2deg(aligned_DoLP_POL(:)), 1), prctile(rad2deg(aligned_DoLP_POL(:)), 99)]));  % [0, 255]
%     bw_AoLP_smooth = imfilter(bw_AoLP, fspecial('average', [3 3]));
%     figure('Name', ['Reflectance + Smoothed AoLP (Grayscale) for ' mainDatasetName]);
%     imshow(reflectance_adj); colormap(jet); colorbar; axis image; hold on;
%     h = imshow(bw_AoLP_smooth, [0 1]); 
%     set(h, 'AlphaData', 0.5);  % Semi-transparent overlay
%     title('AoLP (BW + Brightened + Smoothed)');

end

%% Save fused DoLP image (without markers)
fusedDoLP_fn = fullfile(outputFolder, [mainDatasetName '_fused_DoLP.png']);
saveas(gcf, fusedDoLP_fn);
fprintf('Saved fused DoLP (no markers) to %s\n', fusedDoLP_fn);

%% Selection of points to plot reflectance spectra

% User selection of points
disp('Select points on the fused image. Press Enter when finished.');
% Show fused image for point selection
%hFigPts = figure('Name','Fused DoLP - Select points (click, press Enter when done)');
%imshow(fig); axis image;
title('Click points of interest (press Enter when done)');
hold on;

% User selects points
[x_pts, y_pts] = getpts(); % returns column, row coordinates
x_pts = round(x_pts); y = round(y_pts); % Round coordinates to nearest pixel indices
numPoints = numel(x_pts);

%% Validate points and clamp to image bounds
[H, W, ~] = size(DoLP_POL_rgb);
validIdx = true(1,numPoints);
for k = 1:numPoints
    if x_pts(k) < 1 || x_pts(k) > W || y_pts(k) < 1 || y_pts(k) > H
        warning('Point (%d,%d) out of bounds and will be ignored.', x_pts(k), y_pts(k));
        validIdx(k) = false;
    end
end
x_pts = x_pts(validIdx); y_pts = y_pts(validIdx);
numPoints = numel(x_pts);

% Draw visible red X markers at selected points and save
for k = 1:numPoints
    plot(x_pts(k), y_pts(k), 'x', 'MarkerEdgeColor', 'r', 'MarkerSize', 12, 'LineWidth', 2);
    % Place label next to each pixel
    text(x_pts(k) + 3, y_pts(k) - 3, sprintf('P%d', k), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
end
drawnow;

% Save fused DoLP image with markers
frame = getframe(fusedFig);
img_with_markers = frame.cdata;
withMarkers_fn = fullfile(outputFolder, [mainDatasetName '_fused_DoLP_with_markers.png']);
imwrite(img_with_markers, withMarkers_fn);
fprintf('Saved fused DoLP (with markers) to %s\n', withMarkers_fn);

%% Extract reflectance spectra from selected points

% If there are no selected points, warn and continue
if numPoints == 0
    warning('No valid points selected. Skipping point spectra extraction.');
else  % Extract reflectance spectra

    spectra = zeros(numPoints, size(cube,3));

    % Plot reflectance spectra for selected points
    figure('Name','Reflectance Spectra of Selected Points'); hold on;
    colors = lines(numel(x_pts));

    for i = 1:numPoints
        
        %spectra(i,:) = squeeze(cube(y_pts(i), x_pts(i), :)); %doesn't work as x and y are not positive integers 
        % Convert the point coordinates to valid pixel indices (min/max clamp to the image boundaries to avoid out-of-range access)
        rr = max(1, min(size(cube,1), round(y_pts(i))));
        cc = max(1, min(size(cube,2), round(x_pts(i))));
        spectra(i,:) = squeeze(cube(rr, cc, :)).'; % Extract the spectral vector at the (row rr, column cc) pixel; squeeze removes the singleton dimension, giving an nBandsx1 vector, while the transpose converts it into a 1xnBands row for storage in spectra
   
        legLabels = arrayfun(@(i) sprintf('P%d', i), 1:numPoints, 'UniformOutput', false); % create legend labels for each pixel
        plot(wavelengths, spectra(i,:), 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Point (%d,%d)', x_pts(i), y_pts(i))); hold on;
    end
    axis([400 1000 -0.1 1.1]);
    xlabel('\lambda (nm)'); ylabel('Normalized reflectance (-)');
    title('Reflectance Spectra at Selected Points');
    %legend('show');
    legend(legLabels, 'Location', 'best');
    grid on;

    %% Calculate mean and std reflectance from selected points
    mean_ref = mean(spectra, 1);
    std_ref = std(spectra, 0, 1);
    
    % Plot reflectance spectra with std bounds for selected points
    figure('Name','Mean Reflectance of Selected Points'); hold on;
    plot(wavelengths, mean_ref, 'k-', 'LineWidth', 1.25);
    hold on;
    plot(wavelengths, mean_ref + std_ref, 'r--', 'LineWidth', 1);
    plot(wavelengths, mean_ref - std_ref, 'r--', 'LineWidth', 1);
    axis([400 1000 -0.1 1.1]);
    xlabel('\lambda (nm)'); ylabel('Normalized reflectance (-)');
     legend('Mean','+1 std','-1 std','Location','best');
    title('Mean Reflectance from Selected Points');
    legend('show');
    grid on;

end


%% User selects ROIs for spectral comparison
pairIndex = 1; % allows multiple pairs selection
continueFlag = true;

while continueFlag
        fprintf('\n--- ROI Pair %d ---\n', pairIndex);
    % User selects 2 ROIs at a time for spectral comparison
    reflectances = cell(2,2); % Row1: dataset names, Row2: mean reflectances
    for k = 1:2
        % Show mean image for ROI selection
        figure('Name', sprintf('ROI Pair %d - Draw ROI %d', pairIndex, k));
        imshow(mean(cube,3)); axis image; colormap gray;
        title(sprintf('Draw ROI for dataset %d (press Enter when done)', k));
    
        %h_rect = imrect(); %used in reference code, but not advised by MATLAB
        %pos_rect = round(h_rect.getPosition()); %[x y w h]
        %close(gcf);
    
        h_rect = drawrectangle('Label',sprintf('ROI %d, Pair %d',k, pairIndex),'Color','yellow'); %improved selection instead of imrect
        wait(h_rect); % wait for enter click to finish
        pos_rect = round(h_rect.Position); % [x y w h]
    
    
        % Extract ROI cube
        I_roi = cube(pos_rect(2):(pos_rect(2)+pos_rect(4)), ...
                     pos_rect(1):(pos_rect(1)+pos_rect(3)), :);
    
        % Compute mean and std reflectance in ROI
        [mean_ref, std_ref] = extract_roi_mean(I_roi);
    
        % Store reflectance
        reflectances{1,k} = sprintf('ROI %d (Pair %d)', k, pairIndex);
        reflectances{2,k} = mean_ref(:)'; % 1xN
    
        % Plot mean reflectance with std bounds
        figure('Name',['Reflectance Spectrum - ' reflectances{1,k}]);
        plot(wavelengths, mean_ref, 'b-', 'LineWidth',1.5); hold on;
        plot(wavelengths, mean_ref + std_ref, 'r--', 'LineWidth',1);
        plot(wavelengths, mean_ref - std_ref, 'r--', 'LineWidth',1);
        xlabel('\lambda (nm)'); ylabel('Reflectance (-)');
        title(sprintf('ROI %d, Pair %d', k, pairIndex)); grid on; axis([400 1000 -0.1 1.1]);
    end
    
    % Compute metrics for this ROI pair
    [sid_out, sam_out, sidsam_out, jmsam_out, ns3_out] = compute_spectral_metrics(reflectances{2, 1}, reflectances{2, 2});
    
    % Display Results
    fprintf('\nDatasets: %s and %s\n', reflectances{1, 1}, reflectances{1, 2});
    fprintf('SID: %f\nSAM: %f\nSID-SAM: %f\nJM-SAM: %f\nNS3: %f\n\n', sid_out, sam_out, sidsam_out, jmsam_out, ns3_out);

    % Ask user if they want to add another pair
    answer = questdlg('Select another ROI pair?', 'Continue?', 'Yes','No','No');
    if strcmp(answer, 'No')
        continueFlag = false;
    else
        pairIndex = pairIndex + 1; %repeat whole cycle
    end
end

disp('All PHSI data processed and analysed.');

%% Save all opened images and command window logs

disp(['Saved on: ' string(datetime("now"))]);
disp('=================================');
diary off;

fprintf('Command window log saved to %s\n', logFile);
%save images
if exist('hsiMethod','var')
    save_phsi_results(outputFolder, mainDatasetName, fusionChoice, polMethod, hsiMethod);
else
    save_phsi_results(outputFolder, mainDatasetName, fusionChoice, polMethod);
end

close all