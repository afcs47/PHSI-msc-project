%show_polarization.m/polarimetry.m + dark calibration + DoLP/AoLP(standard) + DoLP/AoLP (fourier series fit) + Methods comparison

clc; clear all; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), 'pol')); % Get auxiliary matlab functions previously created for polarization data analysis


%%  pol_proc()
show = 1; %change to 1 if you want to see the RAW image and the separate polarization images
save = 1; %change to save or not to save only the polarization images
show_aolp = 1; % show only polarimetry (AoLP) image
save_aolp = 1; % save only polarimetry (AoLP) image
rotated_pol_camera = true; % Needed for future fusion with HSI data - change to 'false' for other cases

outputFolder = uigetdir(pwd, 'Select Output Folder'); % Ask where to save results
if outputFolder == 0
    disp('No output folder selected. Exiting.');
    return
end

%outputFolder = 'C:\Users\afili\Downloads\'; % direct path for all output images
% if ~exist(outputFolder, 'dir')
%         mkdir(outputFolder);
% end

% Dark file for calibration
dark_path = 'dark.raw';

%min and max angle (degrees)
angle_min = -90; 
angle_max = 90; % max 180


%change path to file accordingly
%inputFile = "E:\5o ano\Tese\pol_cam\test\images\old\fel.raw"; % 5013504 x 1
%inputFile = "D:\afili\Transferências\Tese - files\pol data\solutions_calib_17jun.raw";

[fileNames, filePath] = uigetfile('*.raw', 'Select one or more RAW files', 'MultiSelect', 'on');
if isequal(fileNames,0)
    disp('No files selected. Exiting.');
    return
end
if ischar(fileNames)
    fileNames = {fileNames}; % Ensure it's a cell array
end


for idx = 1:numel(fileNames)
filename = fileNames{idx};
inputFile = fullfile(filePath, filename);

% Check raw file size
% info = dir(inputFile);
% filesize = info.bytes; %File size: 5013504 bytes
% fprintf('File size: %d bytes\n', filesize);
% filesize/(2048 * 2448) %= bytes per pixel

% Create a subfolder based on the filename variable to save the data
datasetFolder = fullfile(outputFolder, filename);
% Create the folder if it doesn't exist
if ~exist(datasetFolder, 'dir')
    mkdir(datasetFolder);
end
% Use datasetFolder as the new output directory
outputFolder = datasetFolder;


Z = load_image(inputFile, outputFolder , show, save, dark_path);
if rotated_pol_camera; Z = imrotate(Z, 180); end % Only for PHSI system setup

% Perform demosaicing using the RGGB Bayer pattern
J = demosaic(Z, "rggb");
% Display the demosaiced image
figure("Name","Demosaiced image")
imshow(J)
title("Demosaiced image")

[proc_90, proc_45, proc_135, proc_0] = demosaic_polarization_image(Z, inputFile, outputFolder, show, save); % uses  "Z(1:2:end, 1:2:end)" instead of "blockproc(Z, [2 2], @(block) block.data(1,1))" (used in show_polarization.m)

[DoLP, AoLP] = calculate_polarization(proc_90, proc_45, proc_135, proc_0, filename); % DoLP and AoLP from Standard Stokes parameters computation 

visualize_polarization(DoLP, AoLP, inputFile, outputFolder, angle_min, angle_max, show, show_aolp, save, save_aolp, 'Standard');


%% Fourier Series Method for DoLP and AoLP Calculation

% Stack the demosaiced polarization images into a 4D array: [H, W, C, Angle_Index]
polar_images = cat(4, demosaic(uint8(proc_0), 'rggb'), demosaic(uint8(proc_45), 'rggb'), demosaic(uint8(proc_90), 'rggb'),demosaic(uint8(proc_135), 'rggb'));

% Compute 2nd-order
[DoLP_2nd, AoLP_2nd] = calculate_polarization_fourier_2nd(polar_images, filename);
visualize_polarization(DoLP_2nd, AoLP_2nd, inputFile, outputFolder, angle_min, angle_max, show, show_aolp, save, save_aolp, '2nd Order Fourier');

% Compute 4th-order
[DoLP_4th, AoLP_4th] = calculate_polarization_fourier_4th(polar_images, filename);
visualize_polarization(DoLP_4th, AoLP_4th, inputFile, outputFolder, angle_min, angle_max, show, show_aolp, save, save_aolp, '4th Order Fourier');


%% Comparison of Stokes vs Fourier Methods

% Average across RGB channels for comparison
mean_DoLP_stokes = mean(DoLP, 3);
mean_AoLP_stokes = mean(AoLP, 3);

mean_DoLP2 = mean(DoLP_2nd, 3);
mean_AoLP2 = mean(AoLP_2nd, 3);

mean_DoLP4 = mean(DoLP_4th, 3);
mean_AoLP4 = mean(AoLP_4th, 3);

fprintf('\nStandard: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(mean_DoLP_stokes(:)), max(mean_DoLP_stokes(:)), rad2deg(min(mean_AoLP_stokes(:))), rad2deg(max(mean_AoLP_stokes(:))))
fprintf('\n2nd Order Fourier: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(mean_DoLP2(:)), max(mean_DoLP2(:)), rad2deg(min(mean_AoLP2(:))), rad2deg(max(mean_AoLP2(:))))
fprintf('\n4th Order Fourier: DoLP [%.2f %.2f]; AoLP [%.2f %.2f]\n', min(mean_DoLP4(:)), max(mean_DoLP4(:)), rad2deg(min(mean_AoLP4(:))), rad2deg(max(mean_AoLP4(:))))


% DoLP Comparison
compare_methods_DoLP(mean_DoLP_stokes, 'Standard', mean_DoLP2, 'Fourier 2nd', mean_DoLP4, 'Fourier 4th');

% AoLP Comparison
compare_methods_AoLP(mean_AoLP_stokes, 'Standard', mean_AoLP2, 'Fourier 2nd', mean_AoLP4, 'Fourier 4th', angle_min, angle_max)

%% BW plot

% Normalize to grayscale
%bw_DoLP = uint8(255 * mat2gray(mean_DoLP2));  % [0, 255]
%bw_DoLP = uint8(255 * min(max(mean_DoLP2, 0), 1)); %[0, 1] range

% Stretch between 1st and 99th percentiles to brighten the image
low = prctile(mean_DoLP2(:), 1);
high = prctile(mean_DoLP2(:), 99);
bw_DoLP = uint8(255 * mat2gray(mean_DoLP2, [low, high]));

bw_AoLP = uint8(255 * mat2gray(mean_AoLP2));  % [0, 255]

% Show results
fig = figure('Name', 'DoLP & AoLP (Grayscale)');
subplot(1,2,1);
imshow(bw_DoLP);
title('DoLP (2nd Order, BW + Brightened)');

subplot(1,2,2);
imshow(bw_AoLP);
title('AoLP (2nd Order, BW)');

%% Apply smoothing filter
bw_DoLP_smooth = imgaussfilt(bw_DoLP, 2); % sigma = 2 (adjustable)
bw_AoLP_smooth = imgaussfilt(bw_AoLP, 2); %too blurry

kernel = fspecial('average', [2 2]); % [5 5] a bit more blurry
bw_DoLP_smooth = imfilter(bw_DoLP, kernel);
bw_AoLP_smooth = imfilter(bw_AoLP, kernel);


% Show results
fig = figure('Name', 'Smoothed DoLP & AoLP (Grayscale)');
subplot(1,2,1);
imshow(bw_DoLP_smooth);
title('DoLP (2nd Order, BW + Brightened + Smoothed)');

subplot(1,2,2);
imshow(bw_AoLP_smooth);
title('AoLP (2nd Order, BW + Smoothed)');

%% Save results
save_polarization_results(outputFolder, filename, DoLP, AoLP, DoLP_2nd, AoLP_2nd, DoLP_4th, AoLP_4th, bw_DoLP, bw_AoLP, bw_DoLP_smooth, bw_AoLP_smooth);

end
