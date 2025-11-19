%Load all of the files + compute and show DoLP

close all   % Closes all open figure windows 
clear   % Clears all variables from the workspace
clc         % Clears the Command Window

%% Define available datasets
options = {
    'HSI System Only', 'test_datahsi2503_2025-03-25_14-07-36', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2503_2025-03-25_14-07-36\capture\';
    'Polarizer 0', 'test_datahsi2503pol0deg_2025-03-25_14-23-12', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2503pol0deg_2025-03-25_14-23-12\capture\';
    'Polarizer 45', 'test_datahsi2503pol45deg_2025-03-25_14-29-40', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2503pol45deg_2025-03-25_14-29-40\capture\';
    'Polarizer 90', 'test_datahsi2503pol90deg_2025-03-25_14-32-23', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2503pol90deg_2025-03-25_14-32-23\capture\';
    'Polarizer 135', 'test_datahsi2503pol135deg_2025-03-25_14-35-59', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2503pol135deg_2025-03-25_14-35-59\capture\';
    'WG Polarizer 0', 'test_datahsi2603wg0degv2_2025-03-26_10-10-51', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2603wg0degv2_2025-03-26_10-10-51\capture\';
    'WG Polarizer 45', 'test_datahsi2603wg45degv2_2025-03-26_10-15-51', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2603wg45degv2_2025-03-26_10-15-51\capture\';
    'WG Polarizer 90', 'test_datahsi2603wg90deg_2025-03-26_10-21-41', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2603wg90deg_2025-03-26_10-21-41\capture\';
    'WG Polarizer 135', 'test_datahsi2603wg135deg_2025-03-26_10-29-03', 'E:\5o ano\Tese\hsi_cam\test1mar25e26\test_datahsi2603wg135deg_2025-03-26_10-29-03\capture\'};


colors = lines(length(options));
mean_reflectances = struct();
wavelengths = [];

for k = 1:length(options)
    datasetname = options{k, 1};
    fname_selected = options{k, 2};
    fpath = options{k, 3};

    %% Read HS image
    % Read the hyperspectral data and associated information
    [Data, info] = enviread(strcat(fpath,fname_selected,'.raw'),strcat(fpath,fname_selected,'.hdr')); 
    % Read the white reference data
    [White, winfo] = enviread(strcat(fpath,'WHITEREF_',fname_selected,'.raw'),strcat(fpath,'WHITEREF_',fname_selected,'.hdr')); 
    % Read the dark reference data
    [Dark, dinfo] = enviread(strcat(fpath,'DARKREF_',fname_selected,'.raw'),strcat(fpath,'DARKREF_',fname_selected,'.hdr'));
    
    % Extract wavelength information and convert it from string to numeric
    wavelengths = str2num(info.Wavelength(2:end-1));
    
    %% Calibration
    HS_calibrated = zeros(size(Data)); % Initializes a zero matrix for the calibrated data
    white_ref = mean(White, 1); % Averages the white reference across the first dimension
    dark_ref = mean(Dark, 1); % Averages the dark reference across the first dimension
    
    % Calibrates the hyperspectral data by applying the standard formula
    for i = 1:size(Data, 1)
       HS_calibrated(i,:,:) = (Data(i,:,:) - dark_ref(1,:,:))./(white_ref(1,:,:) - dark_ref(1,:,:));
    end

    %% Show fake RGB image
    binning = 8; % Set original binning factor, used for downsampling
    
    % Calculate the approximate indices for the RGB channels
    R = round(333/binning); % channel number (not wavelenght) / binning factor;
    G = round(205/binning); % channel number (not wavelenght) / binning factor;
    B = round(100/binning); % channel number (not wavelenght) / binning factor;
    
    %% Savitzky-Golay Filtering
    hsfiltered = zeros(size(HS_calibrated));  % Initializes the filtered data matrix
    order = 2; % Set order for SG filter
    window = 15; % Set size of window of spectral channels for SG filter 
    
    % Apply the filter to each pixel's spectral data
    for x = 1:length(HS_calibrated(:,1,1))
        for y = 1:length(HS_calibrated(1,:,1))
            hsfiltered(x, y, :) = sgolayfilt(HS_calibrated(x, y, :), order, window);
        end
    end

    %% Selection of region of interest (ROI)
    if k == 1 % First dataset determines the ROI
        validROI = false; % Flag to track if ROI is valid
        while ~validROI
            f = figure('WindowState', 'maximized');
            imshow(hsfiltered(:,:,[R, G, B]),[]) % Displays the filtered image for ROI selection
            title('Select ROI for all datasets');
            h_rect = imrect(); % Allows the user to draw a rectangular ROI
            pos_rect = round(h_rect.getPosition()); % Captures and rounds the ROI position -> [x (leftmost pixel), y (topmost pixel), width, height]
            close(f)

            % Check if the ROI is within bounds
            if pos_rect(2) + pos_rect(4) <= size(hsfiltered,1)... %% Check height limit
                && pos_rect(1) + pos_rect(3) <= size(hsfiltered,2)... % Check width limit
                && pos_rect(2) >= 1 && pos_rect(1) >= 1 % Check non-negative index: top-left corner is not above row 1 and top-left corner is not left of column 1
                
                validROI = true; % Exit loop if ROI is valid
            else
                warndlg('Selected ROI is out of bounds. Please select again.', 'Invalid ROI');
            end
        end
    end

    % Extracts the ROI data from the hs image
    I_roi = hsfiltered(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)),:); % area of image

    %% Calculate mean reflectance from ROI
    sz = size(I_roi); % Gets the size of the ROI data
    
    I = reshape(I_roi, [sz(1) * sz(2), sz(3)]); % Reshapes ROI data into 2D for analysis (pixels x spectral bands)
    
    % Calculates the mean and standard deviation of reflectance across the ROI
    mean_ref = mean(I);
    %std_ref = std(I); 
   
    mean_reflectances.(strrep(datasetname, ' ', '_')) = mean_ref; %add values to respective dataset name, switch spaces to underscore for compatibility
end

%% Compute Stokes parameters (Check dataset existence)
try
    S0_polarizer = mean_reflectances.Polarizer_0 + mean_reflectances.Polarizer_90;
    S1_polarizer = mean_reflectances.Polarizer_0 - mean_reflectances.Polarizer_90;
    S2_polarizer = mean_reflectances.Polarizer_45 - mean_reflectances.Polarizer_135;
    DoLP_polarizer = (sqrt(S1_polarizer.^2 + S2_polarizer.^2) ./ S0_polarizer);
catch
    warning('Polarizer datasets missing!');
    DoLP_polarizer = [];
end

try
    S0_wg = mean_reflectances.WG_Polarizer_0 + mean_reflectances.WG_Polarizer_90;
    S1_wg = mean_reflectances.WG_Polarizer_0 - mean_reflectances.WG_Polarizer_90;
    S2_wg = mean_reflectances.WG_Polarizer_45 - mean_reflectances.WG_Polarizer_135;
    DoLP_wg = (sqrt(S1_wg.^2 + S2_wg.^2) ./ S0_wg);
catch
    warning('Wire Grid Polarizer datasets missing!');
    DoLP_wg = [];
end


%% Compute and Display DoLP max/min
if ~isempty(DoLP_polarizer)
    [maxDOLP_pol, idxMax_pol] = max(DoLP_polarizer);
    [minDOLP_pol, idxMin_pol] = min(DoLP_polarizer);
    fprintf('Polarizer DoLP max: %.4f at %.1f nm\n', maxDOLP_pol, wavelengths(idxMax_pol));
    fprintf('Polarizer DoLP min: %.4f at %.1f nm\n', minDOLP_pol, wavelengths(idxMin_pol));
end

if ~isempty(DoLP_wg)
    [maxDOLP_wg, idxMax_wg] = max(DoLP_wg);
    [minDOLP_wg, idxMin_wg] = min(DoLP_wg);
    fprintf('WG Polarizer DoLP max: %.4f at %.1f nm\n', maxDOLP_wg, wavelengths(idxMax_wg));
    fprintf('WG Polarizer DoLP min: %.4f at %.1f nm\n', minDOLP_wg, wavelengths(idxMin_wg));
end

%% Plot DoLP with Annotations
figure(); hold on;
if ~isempty(DoLP_polarizer)
    plot(wavelengths, DoLP_polarizer, 'r', 'DisplayName', 'DoLP Polarizer');
    scatter([wavelengths(idxMax_pol), wavelengths(idxMin_pol)], [maxDOLP_pol, minDOLP_pol], 'r', 'DisplayName','Max/Min Polarizer')
    text(wavelengths(idxMax_pol), maxDOLP_pol, sprintf('Max: %.4f', maxDOLP_pol), 'Color', 'r');
    text(wavelengths(idxMin_pol), minDOLP_pol, sprintf('Min: %.4f', minDOLP_pol), 'Color', 'r');
end
if ~isempty(DoLP_wg)
    plot(wavelengths, DoLP_wg, 'b', 'DisplayName', 'DoLP WG Polarizer');
    scatter([wavelengths(idxMax_wg), wavelengths(idxMin_wg)], [maxDOLP_wg, minDOLP_wg], 'b', 'DisplayName','Max/Min WG Polarizer')
    text(wavelengths(idxMax_wg), maxDOLP_wg, sprintf('Max: %.4f', maxDOLP_wg), 'Color', 'b');
    text(wavelengths(idxMin_wg), minDOLP_wg, sprintf('Min: %.4f', minDOLP_wg), 'Color', 'b');
end
xlabel('\lambda (nm)'); ylabel('DoLP');
title('Degree of Linear Polarization (DoLP)');
legend show;
grid on; grid minor;
hold off;

%% Plot DoLP with Normalized Axis
figure(); hold on;

if ~isempty(DoLP_polarizer)
    plot(wavelengths, DoLP_polarizer, 'r', 'DisplayName', 'DoLP Polarizer (norm)');
end

if ~isempty(DoLP_wg)
    plot(wavelengths, DoLP_wg, 'b', 'DisplayName', 'DoLP WG Polarizer (norm)');
end

xlabel('\lambda (nm)'); ylabel('Normalized DoLP');
title('Normalized Degree of Linear Polarization (DoLP)');
legend show;
grid on; grid minor;
ylim([0 1]); % force normalization scale
hold off;


%% Compute AoLP
try
    AoLP_polarizer = 0.5 * atan2(S2_polarizer, S1_polarizer); % radians
catch
    warning('Cannot compute AoLP for Polarizer (missing data)');
    AoLP_polarizer = [];
end

try
    AoLP_wg = 0.5 * atan2(S2_wg, S1_wg); % radians
catch
    warning('Cannot compute AoLP for WG Polarizer (missing data)');
    AoLP_wg = [];
end

%% Convert AoLP to degrees for readability
AoLP_polarizer_deg = rad2deg(AoLP_polarizer);
AoLP_wg_deg = rad2deg(AoLP_wg);

%% Plot AoLP
figure(); hold on;
if ~isempty(AoLP_polarizer)
    plot(wavelengths, AoLP_polarizer_deg, 'r', 'DisplayName', 'AoLP Polarizer');
end
if ~isempty(AoLP_wg)
    plot(wavelengths, AoLP_wg_deg, 'b', 'DisplayName', 'AoLP WG Polarizer');
end
xlabel('\lambda (nm)'); ylabel('AoLP (degrees)');
title('Angle of Linear Polarization (AoLP)');
legend show;
ylim([-90,90])
grid on; grid minor;
hold off;

%% Compute and Display DoLP max/min
if ~isempty(AoLP_polarizer)
    [maxAOLP_pol, idxMax_pol] = max(AoLP_polarizer_deg);
    [minAOLP_pol, idxMin_pol] = min(AoLP_polarizer_deg);
    fprintf('Polarizer AoLP max: %.4f deg at %.1f nm\n', maxAOLP_pol, wavelengths(idxMax_pol));
    fprintf('Polarizer AoLP min: %.4f deg at %.1f nm\n', minAOLP_pol, wavelengths(idxMin_pol));
end

if ~isempty(AoLP_wg)
    [maxAOLP_wg, idxMax_wg] = max(AoLP_wg_deg);
    [minAOLP_wg, idxMin_wg] = min(AoLP_wg_deg);
    fprintf('WG Polarizer AoLP max: %.4f deg at %.1f nm\n', maxAOLP_wg, wavelengths(idxMax_wg));
    fprintf('WG Polarizer AoLP min: %.4f deg at %.1f nm\n', minAOLP_wg, wavelengths(idxMin_wg));
end