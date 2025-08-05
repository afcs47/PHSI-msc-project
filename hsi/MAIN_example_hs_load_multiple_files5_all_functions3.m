%All functions put together + DoLP spatial distribution + AoLP

close all   % Closes all open figure windows
clear   % Clears all variables from the workspace
clc         % Clears the Command Window

%% Define available datasets
options = {...
    'HSI System Only', 'test_datahsi2503_2025-03-25_14-07-36', 'test_datahsi2503_2025-03-25_14-07-36\capture\';
    'Polarizer 0', 'test_datahsi2503pol0deg_2025-03-25_14-23-12', 'test_datahsi2503pol0deg_2025-03-25_14-23-12\capture\';
    'Polarizer 45', 'test_datahsi2503pol45deg_2025-03-25_14-29-40', 'test_datahsi2503pol45deg_2025-03-25_14-29-40\capture\';
    'Polarizer 90', 'test_datahsi2503pol90deg_2025-03-25_14-32-23', 'test_datahsi2503pol90deg_2025-03-25_14-32-23\capture\';
    'Polarizer 135', 'test_datahsi2503pol135deg_2025-03-25_14-35-59', 'test_datahsi2503pol135deg_2025-03-25_14-35-59\capture\';
    'WG Polarizer 0', 'test_datahsi2603wg0degv2_2025-03-26_10-10-51', 'test_datahsi2603wg0degv2_2025-03-26_10-10-51\capture\';
    'WG Polarizer 45', 'test_datahsi2603wg45degv2_2025-03-26_10-15-51', 'test_datahsi2603wg45degv2_2025-03-26_10-15-51\capture\';
    'WG Polarizer 90', 'test_datahsi2603wg90deg_2025-03-26_10-21-41', 'test_datahsi2603wg90deg_2025-03-26_10-21-41\capture\';
    'WG Polarizer 135', 'test_datahsi2603wg135deg_2025-03-26_10-29-03', 'test_datahsi2603wg135deg_2025-03-26_10-29-03\capture\' };

pos_rect = zeros(1, 4);

while true 
    %% Ask user to select an option
    %choice_initial = menu("Select an option:", ["Spectra plotting" "DoLP and AoLP computation" "Spectral metrics calculation"]);
    choice_initial = menu("Select an option:", "Spectra plotting", "DoLP and AoLP computation", "Spectral metrics calculation", 'Exit');
    
    if (choice_initial == 1)
        % Ask user to select multiple datasets
        [selectedIndexes, tf] = listdlg('PromptString','Select datasets:', 'SelectionMode','multiple', 'ListString', options(:,1));
        if ~tf
            continue; % Return to main menu
        end
        % Ask user for plotting settings
        choice_std_dev = questdlg('Consider standard deviation?', 'Next Action', 'Yes', 'No', 'No');
        if strcmp(choice_std_dev, 'Yes')
                choice_fill = questdlg('Fill in between min and max?', 'Next Action', 'Yes', 'No', 'No');
        else 
            choice_fill = 'No';
        end
        
        colors = choose_colors(selectedIndexes);
    
        figure(); hold on; % Initialize figure for comparison
        title('Comparison of Datasets');
        axis([400 1000 0 1.2]) % Sets the axis limits
        xlabel('\lambda (nm)') % X-axis label (wavelengths)
        ylabel('Normalized reflectance (-)') % Y-axis label
        grid on;
        grid minor;
    
        for k = 1:length(selectedIndexes)
            selectedIndex = selectedIndexes(k);
            datasetname = options{selectedIndex, 1};
            fname_selected = options{selectedIndex, 2};
            fpath = options{selectedIndex, 3};
    
            [wavelengths, pos_rect, mean_ref, std_ref] = process_data(fpath, fname_selected, datasetname, false, k, pos_rect);
    
            % Plot mean reflectance
            plot_mean_reflectance(datasetname, wavelengths, mean_ref, std_ref, choice_std_dev, choice_fill, colors, k)
        end
        legend show; % Show legend for dataset labels
        hold off;


    elseif (choice_initial==2)
        choice_2 = menu("Select an option:", "DoLP/AoLP as a spectra function", "DoLP/AoLP as a spatial function", "Cancel");
            if (choice_2==1) % as a spectra function
                mean_reflectances = struct();
                for k = 1:length(options)
                    datasetname = options{k, 1};
                    fname_selected = options{k, 2};
                    fpath = options{k, 3};
        
                    [wavelengths, pos_rect, mean_ref, std_ref] = process_data(fpath, fname_selected, datasetname, false, k, pos_rect);
                    mean_reflectances.(strrep(datasetname, ' ', '_')) = mean_ref;
                end

                % Compute Stokes parameters and DoLP
                DoLP_polarizer = compute_dolp(mean_reflectances, 'Polarizer');
                [maxDoLP_pol, minDoLP_pol, idxMax_pol, idxMin_pol] = compute_extremes(DoLP_polarizer, 'Polarizer', wavelengths);
                
                DoLP_wg = compute_dolp(mean_reflectances, 'WG_Polarizer');
                [maxDoLP_wg, minDoLP_wg, idxMax_wg, idxMin_wg] = compute_extremes(DoLP_wg, 'WG_Polarizer', wavelengths);
                
                % Plot DoLP with annotations
                plot_polarization_parameter(wavelengths, DoLP_polarizer, DoLP_wg, idxMax_pol, idxMin_pol, maxDoLP_pol, minDoLP_pol, idxMax_wg, idxMin_wg, maxDoLP_wg, minDoLP_wg, 'DoLP', 'Degree of Linear Polarization (DoLP)');


                % Compute AoLP for each wavelength
                AoLP_polarizer = compute_aolp(mean_reflectances, 'Polarizer');
                [maxAoLP_pol, minAoLP_pol, idxMax_pol, idxMin_pol] = compute_extremes(AoLP_polarizer, 'Polarizer', wavelengths);
                
                AoLP_wg = compute_aolp(mean_reflectances, 'WG_Polarizer');
                [maxAoLP_wg, minAoLP_wg, idxMax_wg, idxMin_wg] = compute_extremes(AoLP_wg, 'WG_Polarizer', wavelengths);
                
                % Plot AoLP with annotations
                plot_polarization_parameter(wavelengths, AoLP_polarizer, AoLP_wg, idxMax_pol, idxMin_pol, maxAoLP_pol, minAoLP_pol, idxMax_wg, idxMin_wg, maxAoLP_wg, minAoLP_wg, 'AoLP (degrees)', 'Angle of Linear Polarization (AoLP)');


            elseif (choice_2==2) % as a spatial function
                spatial_reflectances = struct();
                for k = 1:length(options)
                    datasetname = options{k, 1};
                    fname_selected = options{k, 2};
                    fpath = options{k, 3};

                    % Extract the calibrated hyperspectral data (HS_calibrated) and region of interest (ROI)
                    [Data, White, Dark, wavelengths] = read_data(fpath, fname_selected); % Read HS image
                    HS_calibrated = apply_calibration(Data, White, Dark); % Calibration
                    [R, G, B] = fake_rgb(false, HS_calibrated, datasetname); % Fake RGB image
                    hsfiltered = apply_sg_filter(HS_calibrated); % Savitzky-Golay Filtering
                    
                    % Compute mean reflectance over the wavelengths
                    mean_ref = mean(hsfiltered, 3); % x, y, lambda

                    % Store results
                    spatial_reflectances.(strrep(datasetname, ' ', '_')) = mean_ref;
 
                end
       
                % Compute DoLP for each pixel
                DoLP_polarizer = compute_dolp(spatial_reflectances, 'Polarizer');
                DoLP_wg = compute_dolp(spatial_reflectances, 'WG_Polarizer');
        
                plot_polarization_spatially(DoLP_polarizer, DoLP_wg, 'turbo', 'Spatial Distribution of DoLP');


                % Compute AoLP for each pixel
                AoLP_polarizer = compute_aolp(spatial_reflectances, 'Polarizer');
                AoLP_wg = compute_aolp(spatial_reflectances, 'WG_Polarizer');

                plot_polarization_spatially(AoLP_polarizer, AoLP_wg, 'turbo', 'Spatial Distribution of AoLP');
            else
               continue;
            end

    elseif (choice_initial==3) 
    % Calculate some spectral measures (additional (Add on) hyperspectral toolbox for image processing toolbox necessary)
    % https://www.mathworks.com/help/images/hyperspectral-image-processing.html
    
    % Selects two sets of mean reflectance from user-defined ROIs
        
        cancel_selection = false;  % Flag to track cancellation
        reflectances = cell(2, 2); % Row 1 - dataset names; Row 2 - reflectance values
        for k=1:2
            % Ask user to select a dataset
            [selectedIndex, tf] = listdlg('PromptString','Select dataset:', 'SelectionMode','single', 'ListString', options(:,1));
            if ~tf
                cancel_selection = true; % User canceled
                break; % Exit dataset selection loop
            end

            % Assign selected dataset variables
            datasetname = options{selectedIndex, 1};
            fname_selected = options{selectedIndex, 2};
            fpath = options{selectedIndex, 3};
    
            [wavelengths, ~, mean_ref, std_ref] = process_data(fpath, fname_selected, datasetname, true, 1, pos_rect);
            
            %size(mean_ref)
            %reflectances(:, k) = mean_ref';  

            % Store dataset name
            reflectances{1, k} = datasetname;  
            % Store reflectance values
            reflectances{2, k} = mean_ref';  % Transpose mean_ref since it's 1×N


            % Plot mean reflectance
            figure()
            plot(wavelengths, mean_ref,'b-') % Plots the mean reflectance in blue
            hold on
            plot(wavelengths, mean_ref + std_ref, 'r-') % Plots the upper bound (mean + std)
            plot(wavelengths, mean_ref - std_ref, 'r-') % Plots the lower bound (mean - std)
            axis([400 1000 0 1.2]) % Sets the axis limits
            xlabel('\lambda (nm)') % X-axis label (wavelengths)
            ylabel('Normalized reflectance (-)') % Y-axis label
            title(datasetname)
        end

        % If the user canceled, returns to the main menu
        if cancel_selection
            continue;
        end

        [sid_out, sam_out, sidsam_out, jmsam_out, ns3_out] = compute_spectral_metrics(reflectances{2, 1}, reflectances{2, 2});

        % Display Results
        fprintf('\nDatasets: %s and %s\n', reflectances{1, 1}, reflectances{1, 2});
        fprintf('SID: %f\nSAM: %f\nSID-SAM: %f\nJM-SAM: %f\nNS3: %f\n', sid_out, sam_out, sidsam_out, jmsam_out, ns3_out);

    else
       return;
    end
end




%% Functions:

function [Data, White, Dark, wavelengths] = read_data(fpath, fname)
    [Data, info] = enviread(strcat(fpath,fname,'.raw'), strcat(fpath,fname,'.hdr')); % Read the hyperspectral data and associated information
    [White, ~] = enviread(strcat(fpath,'WHITEREF_',fname,'.raw'), strcat(fpath,'WHITEREF_',fname,'.hdr')); % Read the white reference data
    [Dark, ~] = enviread(strcat(fpath,'DARKREF_',fname,'.raw'), strcat(fpath,'DARKREF_',fname,'.hdr')); % Read the dark reference data
    wavelengths = str2num(info.Wavelength(2:end-1)); % Extract wavelength information and convert it from string to numeric
end


function colors = choose_colors(n)
    if (length(n)<= 6)
       colors = lines(length(n)); % Generate distinct colors for plots
    else 
       colors = turbo(length(n));
    end
end


function HS_calibrated = apply_calibration(Data, White, Dark)
    HS_calibrated = zeros(size(Data)); % Initializes a zero matrix for the calibrated data
    white_ref = mean(White, 1); % Averages the white reference across the first dimension
    dark_ref = mean(Dark, 1); % Averages the dark reference across the first dimension
    
    % Calibrates the hyperspectral data by applying the standard formula
    for i = 1:size(Data, 1)
       HS_calibrated(i,:,:) = (Data(i,:,:) - dark_ref(1,:,:))./(white_ref(1,:,:) - dark_ref(1,:,:));
    end
end


function [R, G, B] = fake_rgb(showRGB, HS_calibrated, datasetname)
    binning = 8; % Set original binning factor, used for downsampling
    
    % Calculate the approximate indices for the RGB channels
    R = round(333/binning); % channel number (not wavelenght) / binning factor;
    G = round(205/binning); % channel number (not wavelenght) / binning factor;
    B = round(100/binning); % channel number (not wavelenght) / binning factor;
    
    if showRGB
        % Display the calibrated data as a pseudo-RGB image
        figure()
        imshow(HS_calibrated(:,:,[R, G, B]),[])
        title(datasetname)
    end 
end


function hsfiltered = apply_sg_filter(HS_calibrated)
    hsfiltered = zeros(size(HS_calibrated));  % Initializes the filtered data matrix
    order = 2; % Set order for SG filter
    window = 15; % Set size of window of spectral channels for SG filter 
    
    % Apply the filter to each pixel's spectral data
    for x = 1:length(HS_calibrated(:,1,1))
        for y = 1:length(HS_calibrated(1,:,1))
            hsfiltered(x, y, :) = sgolayfilt(HS_calibrated(x, y, :), order, window);
        end
    end
end


function [validROI, pos_rect] = select_roi(hsfiltered,R,G,B)
    f = figure('WindowState', 'maximized');
    imshow(hsfiltered(:,:,[R, G, B]),[]) % Displays the image for ROI selection
    title('Select ROI:')
    h_rect = imrect(); % Allows the user to select the ROI
    pos_rect = round(h_rect.getPosition()); % Gets and rounds the ROI position
    close(f) % Closes the figure

    % Check if the ROI is within bounds
    if pos_rect(2) + pos_rect(4) <= size(hsfiltered,1)... %% Check height limit
        && pos_rect(1) + pos_rect(3) <= size(hsfiltered,2)... % Check width limit
        && pos_rect(2) >= 1 && pos_rect(1) >= 1 % Check non-negative index: top-left corner is not above row 1 and top-left corner is not left of column 1
        
        validROI = true; % Exit loop if ROI is valid
    else
        warndlg('Selected ROI is out of bounds. Please select again.', 'Invalid ROI');
        validROI = false;
    end
end


function [mean_ref_fc, std_ref] = extract_roi_mean(I_roi)
    sz = size(I_roi); 
    I = reshape(I_roi, [sz(1) * sz(2), sz(3)]); % Reshapes data for processing
    mean_ref_fc = mean(I); % Computes the mean reflectance
    std_ref = std(I); % Computes the standard deviation
end


function [wavelengths, pos_rect, mean_ref, std_ref] = process_data(fpath, fname_selected, datasetname, showRGB, k, pos_rect)
    % Read HS image
    [Data, White, Dark, wavelengths] = read_data(fpath, fname_selected);
    
    % Calibration
    HS_calibrated = apply_calibration(Data, White, Dark);
    
    % Show fake RGB image
    [R, G, B] = fake_rgb(showRGB, HS_calibrated, datasetname);
    
    % Savitzky-Golay Filtering
    hsfiltered = apply_sg_filter(HS_calibrated);

    % Selection of region of interest (ROI)
    if k == 1 % First dataset determines the ROI
        validROI = false; % Flag to track if ROI is valid
        while ~validROI
            [validROI, pos_rect] = select_roi(hsfiltered,R,G,B);
        end
    end
    I_roi = hsfiltered(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)),:); % area of image

    % Extract the ROI data from the hs image and calculate mean reflectance from it
    [mean_ref, std_ref] = extract_roi_mean(I_roi);
end


function plot_mean_reflectance(datasetname, wavelengths, mean_ref, std_ref, choice_std_dev, choice_fill, colors, k)
    plot(wavelengths, mean_ref,'Color', colors(k,:), 'DisplayName', datasetname); % Plots the mean reflectance for the given dataset      
    if strcmp(choice_std_dev, 'Yes')
        if strcmp(choice_fill, 'No')
            plot(wavelengths, mean_ref + std_ref, '-', 'Color', colors(k,:), 'DisplayName', strcat(datasetname, ' Max')); % Plots the upper bound (mean + std)
            plot(wavelengths, mean_ref - std_ref, '-', 'Color', colors(k,:), 'DisplayName', strcat(datasetname, ' Min')); % Plots the lower bound (mean - std)
        elseif strcmp(choice_fill, 'Yes')
            % Plot mean reflectance with shaded standard deviation
            x_vals = [wavelengths, fliplr(wavelengths)]; % X values for fill (wavelengths mirrored)
            y_vals = [mean_ref - std_ref, fliplr(mean_ref + std_ref)]; % Y values for fill
            
            % Plot shaded region for standard deviation
            fill(x_vals, y_vals, colors(k,:), 'FaceAlpha', 0.2, 'EdgeColor', colors(k,:), 'DisplayName', strcat(datasetname, ' Std Dev')); 
        end
    end 
end


function DoLP = compute_dolp(mean_reflectances, type)
    % Compute Stokes parameters
    S0 = mean_reflectances.([type '_0']) + mean_reflectances.([type '_90']); % S0 = L0 + L90
    S1 = mean_reflectances.([type '_0']) - mean_reflectances.([type '_90']); % S1 = L0 - L90
    S2 = mean_reflectances.([type '_45']) - mean_reflectances.([type '_135']); % S2 = L45 - L135
    
    % Compute Degree of Linear Polarization (DoLP)
    DoLP = sqrt(S1.^2 + S2.^2) ./ S0; 
end


function [max_x, min_x, idxMax, idxMin] = compute_extremes(x, type, wavelengths) 
    % Find max and min DoLP
    [max_x, idxMax] = max(x);
    [min_x, idxMin] = min(x);
    
    % Extract variable name before '_'
    varName = regexp(inputname(1), '^[^_]+', 'match', 'once'); % ^ = start of the string; [^_] = Any character except _ ; + = One or more occurrences
    
    % Display results
    fprintf('\n%s %s max: %.4f at %.1f nm\n', type, varName, max_x, wavelengths(idxMax));
    fprintf('%s %s min: %.4f at %.1f nm\n', type, varName, min_x, wavelengths(idxMin));
end



function [sid_out, sam_out, sidsam_out, jmsam_out, ns3_out] = compute_spectral_metrics(reflectance1,reflectance2)

    % Calculate some spectral measures (additional (Add on) hyperspectral toolbox for image processing toolbox necessary)
    % https://www.mathworks.com/help/images/hyperspectral-image-processing.html
    
    % Calculates various spectral similarity metrics
    sid_out = sid(reflectance1,reflectance2); % Spectral Information Divergence
    sam_out = sam(reflectance1,reflectance2); % Spectral Angle Mapper
    sidsam_out = sidsam(reflectance1,reflectance2); % Combined SID and SAM metric
    jmsam_out = jmsam(reflectance1,reflectance2); % Jeffries Matusita-Spectral Angle Mapper 
    ns3_out = ns3(reflectance1,reflectance2); % Normalized Spectral Similarity Score
end

function AoLP = compute_aolp(mean_reflectances, type)
    % Compute Stokes parameters
    S1 = mean_reflectances.([type '_0']) - mean_reflectances.([type '_90']); % S1 = L0 - L90
    S2 = mean_reflectances.([type '_45']) - mean_reflectances.([type '_135']); % S2 = L45 - L135

    % Compute Angle of Linear Polarization (AoLP)
    AoLP = 0.5 * atan2(S2, S1); % atan2 for correct quadrant handling

    % Convert radians to degrees
    AoLP = rad2deg(AoLP);
end


function plot_polarization_parameter(wavelengths, data1, data2, idxMax1, idxMin1, max1, min1, idxMax2, idxMin2, max2, min2, yLabel, titleText)
    figure(); hold on;
    
    % Plot first dataset - polarizer
    plot(wavelengths, data1, 'r', 'DisplayName','Polarizer');
    scatter([wavelengths(idxMax1), wavelengths(idxMin1)], [max1, min1], 'r', 'DisplayName', 'Max/Min Polarizer');
    text(wavelengths(idxMax1), max1, sprintf('Max: %.4f', max1), 'Color', 'r');
    text(wavelengths(idxMin1), min1, sprintf('Min: %.4f', min1), 'Color', 'r');

    % Plot second dataset - wg
    plot(wavelengths, data2, 'b', 'DisplayName','WG Polarizer');
    scatter([wavelengths(idxMax2), wavelengths(idxMin2)], [max2, min2], 'b', 'DisplayName', 'Max/Min WG Polarizer');
    text(wavelengths(idxMax2), max2, sprintf('Max: %.4f', max2), 'Color', 'b');
    text(wavelengths(idxMin2), min2, sprintf('Min: %.4f', min2), 'Color', 'b');

    % Labels, title, and formatting
    xlabel('\lambda (nm)');
    ylabel(yLabel);
    title(titleText);
    legend show;
    grid on; grid minor;
    
    hold off;
end


function plot_polarization_spatially(data1, data2, colormap_name, titleText)   
    % Compute common color scale limits
    Lim = [min([data1(:); data2(:)]), max([data1(:); data2(:)])];
    %Lim = [0, 1]; %DoLP
    %Lim = [-90, 90]; %AoLP

    figure();
    % Create tiled layout for better visualization
    tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % First image - polarizer
    nexttile;
    imagesc(data1);
    colorbar;
    colormap(colormap_name);
    title([titleText, ' - Polarizer']); 
    xlabel('X Position'); 
    ylabel('Y Position');
    axis image; % Ensure square pixels
    clim(Lim); % Set a common color scale

    % Second image - wg polarizer
    nexttile;
    imagesc(data2);
    colorbar;
    colormap(colormap_name);
    title([titleText, ' - WG Polarizer']);
    xlabel('X Position'); 
    ylabel('Y Position');
    axis image; % Ensure square pixels
    clim(Lim); % Set a common color scale
end


