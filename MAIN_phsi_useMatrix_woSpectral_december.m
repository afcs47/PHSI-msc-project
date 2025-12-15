%% MAIN SCRIPT FOR HSI + POLARIZATION FUSION and Overlay Analysis - use matrix previously computed for calibration and skip spectral analysis (ie just fusion plots)

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
        %tform = calibrate_hsi_polData_from_checkerboard(hsi_rgb_u8, pol_rgb, [selectedType '_calib_' selectedDay], size(DoLP_HSI), size(DoLP_POL), check_calibration);
        % Build expected filename in same folder as this script
        scriptDir = fileparts(mfilename('fullpath'));
        hFile = fullfile(scriptDir, sprintf('%s_calib_%s_homography_tform.mat', selectedType, selectedDay));
        
        if exist(hFile, 'file')
            fprintf('Homography file found. Loading existing matrix: %s\n', hFile);
            S = load(hFile);
            
            if isfield(S, 'tform')
                tform = S.tform;
                fprintf('Loaded saved homography.\n');
            else
                warning('File exists but does not contain variable "tform". Recomputing homography...');
                tform = calibrate_hsi_polData_from_checkerboard(hsi_rgb_u8, pol_rgb, [selectedType '_calib_' selectedDay], size(DoLP_HSI), size(DoLP_POL), check_calibration);
                save(hFile, 'tform');
                fprintf('Computed and saved new homography to %s\n', hFile);
            end
        end
    case 'HSI without WG + POL'
        % Only reflectance fusion possible
        [reflectance_HSI, ~] = load_hsi_analysis('Reflectance', selectedType, [selectedType '_' woWGCalibData.Day{1} '_woWG'], 'Spatial_reflectances'); %selectedDay of sampled data = day of calibration data
        [DoLP_POL, AoLP_POL] = load_pol_analysis(polMethod, [selectedType '_woWG_' woWGCalibData.Day{1}]);
        [cube, wavelengths] = load_hsi_analysis('Reflectance', selectedType, [selectedType '_' woWGCalibData.Day{1} '_woWG'], 'Spectral_reflectances');
        reflectance_HSI = reflectance_HSI.sample_01;
        cube = cube.sample_01;
        %% Calibration: Get homography from checkerboard
        %tform = calibrate_hsi_polData_from_checkerboard(hsi_rgb_u8, pol_rgb, [selectedType '_woWG_calib_' selectedDay], size(reflectance_HSI), size(DoLP_POL), check_calibration);
        % Build expected filename in same folder as this script
        scriptDir = fileparts(mfilename('fullpath'));
        hFile = fullfile(scriptDir, sprintf('%s_woWG_calib_%s_homography_tform.mat', selectedType, selectedDay));
        
        if exist(hFile, 'file')
            fprintf('Homography file found. Loading existing matrix: %s\n', hFile);
            S = load(hFile);
            
            if isfield(S, 'tform')
                tform = S.tform;
                fprintf('Loaded saved homography.\n');
            else
                warning('File exists but does not contain variable "tform". Recomputing homography...');
                tform = calibrate_hsi_polData_from_checkerboard(hsi_rgb_u8, pol_rgb, [selectedType '_woWG_calib_' selectedDay], size(DoLP_HSI), size(DoLP_POL), check_calibration);
                save(hFile, 'tform');
                fprintf('Computed and saved new homography to %s\n', hFile);
            end
        end
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


end

%% Save fused DoLP image (without markers)
fusedDoLP_fn = fullfile(outputFolder, [mainDatasetName '_fused_DoLP.png']);
saveas(gcf, fusedDoLP_fn);
fprintf('Saved fused DoLP (no markers) to %s\n', fusedDoLP_fn);

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
close all;