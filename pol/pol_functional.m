function Z = load_image(show, output_path, input_file, save)
    [~, baseFilename, ~] = fileparts(input_file);

    % load image with given dimensions
    row=2048;  col=2448;
    fin=fopen(input_file,'r');
    I=fread(fin,row*col,'uint8=>uint8'); 
    Z=reshape(I,col,row);
    Z=Z';
        if show == true 
            figure('Visible','on')
            imshow(Z);
            ax = gca;
        else
            figure(Visible="off")
            imshow(Z);
            ax = gca;
            % title("Input RAW image")
        end
        if save == 1 
            exportgraphics(ax, fullfile(output_path, [baseFilename '_raw.png']), 'Resolution', 400);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DoLP, AoLP] = calculate_polarization(Z)
    % Extract four polarization orientations from image
    pol_90 = Z(1:2:end, 1:2:end);  % Top-left pixels
    pol_45 = Z(1:2:end, 2:2:end);  % Top-right pixels
    pol_135 = Z(2:2:end, 1:2:end); % Bottom-left pixels
    pol_0 = Z(2:2:end, 2:2:end);   % Bottom-right pixels
    
    % Demosaic each polarization separately
    color_90 = demosaic(uint8(pol_90), 'rggb');
    color_45 = demosaic(uint8(pol_45), 'rggb');
    color_135 = demosaic(uint8(pol_135), 'rggb');
    color_0 = demosaic(uint8(pol_0), 'rggb');
    
    [height, width, channels] = size(color_0);
    
    % init DoLP and AoLP
    DoLP = zeros(height, width, channels);
    AoLP = zeros(height, width, channels);
    
    % Process each color channel
    for c = 1:channels
        % intensities for each polarization angle (need to be type double to
        % work in sqrt function)
        I_0 = double(color_0(:,:,c));
        I_90 = double(color_90(:,:,c));
        I_45 = double(color_45(:,:,c));
        I_135 = double(color_135(:,:,c));
        
        % Calculate Stokes parameters
        S0 = I_0 + I_90;  % Total image intensity
        S1 = I_0 - I_90;  % diff between horizontal and vertical polarization
        S2 = I_45 - I_135; % diff between 45° and 135° polarization
        
        % Calculate DoLP (Degree of Linear Polarization)
        DoLP(:,:,c) = sqrt(S1.^2 + S2.^2) ./ (S0 + eps);
        
        % Calculate AoLP (Angle of Linear Polarization)
        AoLP(:,:,c) = 0.5 * atan2(S2, S1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_polarization(DoLP, AoLP, output_path, input_file, save)
    % Create HSV image where:
    % - Hue represents AoLP (scaled to account for π periodicity)
    % - Saturation set to 1
    % - Value is set to 1
    [~, baseFilename, ~] = fileparts(input_file);
    [height, width, channels] = size(DoLP);
    
    % Figure for mean DoLP across all channels
    dolp_fig = figure('Name', 'Mean DoLP', 'Visible', 'on');
    meanDoLP = mean(DoLP, 3);
    imshow(meanDoLP);
    colormap(jet);
    colorbar;
    % title('DoLP Across Channels','FontSize',14);
    
    % Figure for mean AoLP
    aolp_fig = figure('Name', 'Mean AoLP', 'Visible', 'on');
    
    % Create axes for the AoLP image that leaves space for the colorwheel
    aolp_ax = axes('Position', [0.05, 0.05, 0.75, 0.9]);
    
    % Calculate mean AoLP
    x_total = zeros(height, width);
    y_total = zeros(height, width);
    
    % For each channel, add its weighted contribution
    for c = 1:channels
        % Convert to Cartesian coordinates to handle the circular nature of angles
        x_component = cos(2 * AoLP(:,:,c));
        y_component = sin(2 * AoLP(:,:,c));
        
        % Add to accumulated components
        x_total = x_total + x_component;
        y_total = y_total + y_component;
    end
    
    % Convert back to angle
    meanAoLP = 0.5 * atan2(y_total, x_total);
    meanAoLP = mod(meanAoLP, pi);
    
    % Create HSV image for mean AoLP
    hsv_img = zeros(height, width, 3);
    hsv_img(:,:,1) = meanAoLP / pi; % Hue from AoLP
    hsv_img(:,:,2) = ones(size(meanAoLP)); % Full saturation
    hsv_img(:,:,3) = ones(size(meanAoLP)); % Value = 1
    
    % Convert to RGB and display
    rgb_img = hsv2rgb(hsv_img);
    imshow(rgb_img, 'Parent', aolp_ax);
    % title(aolp_ax, 'AoLP across channels');   
    
    % Add colorwheel to the AoLP figure with specific position
    polarization_colorwheel(aolp_fig, [0.82, 0.3, 0.15, 0.15]);
    
    % save the images
    if save == 1 %change to false to stop saving
        exportgraphics(dolp_fig, fullfile(output_path, [baseFilename '_dolp.png']), 'Resolution', 400);
        exportgraphics(aolp_fig, fullfile(output_path, [baseFilename '_aolp.png']), 'Resolution', 400);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function polarization_colorwheel(fig, position)
    % Add a color wheel to show the mapping between color and polarization angle
    % Create a separate axes for the color wheel in the specified figure
    if nargin < 2
        position = [0.75, 0.1, 0.2, 0.2]; % Default position if not specified
    end
    
    % Use the parent figure and create a new axes with specified position
    axes('Parent', fig, 'Position', position);
    
    % Create color wheel
    [x, y] = meshgrid(linspace(-1, 1, 100), linspace(-1, 1, 100));
    r = sqrt(x.^2 + y.^2);
    
    % Calculate angle in range [0, π]
    % This accounts for the π periodicity of polarization
    theta = mod(atan2(y, x), pi);
    
    % Convert to HSV
    hue = theta / pi; % Map [0, π] to [0, 1]
    sat = r;
    val = ones(size(r));
    
    % Mask out points outside the unit circle
    mask = r > 1;
    hue(mask) = 0;
    sat(mask) = 0;
    val(mask) = 0;
    
    % Create HSV image
    hsv_img = zeros(100, 100, 3);
    hsv_img(:,:,1) = hue;
    hsv_img(:,:,2) = sat;
    hsv_img(:,:,3) = val;
    
    % Convert to RGB and display
    rgb_wheel = hsv2rgb(hsv_img);
    imshow(rgb_wheel);
    axis off;
    title('AoLP Reference', 'FontSize', 8);
    
    % Add angle labels around the wheel
    hold on;
    angles = [0, pi/4, pi/2, 3*pi/4];
    labels = {'0°', '45°', '90°', '135°'};
    for i = 1:length(angles)
        angle = angles(i);
        x_pos = 0.8 * cos(angle);
        y_pos = 0.8 * sin(angle);
        text(x_pos*50+50, y_pos*50+50, labels{i}, 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 8);
    end
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function demosaic_polarization_image(Z,output_path, input_file, save, show)
[~, baseFilename, ~] = fileparts(input_file);

% Extract polarization images directly
proc1 = Z(1:2:end, 1:2:end);     % 90° (top-left pixels)
proc2 = Z(1:2:end, 2:2:end);     % 45° (top-right pixels)
proc3 = Z(2:2:end, 1:2:end);     % 135° (bottom-left pixels)
proc4 = Z(2:2:end, 2:2:end);     % 0° (bottom-right pixels)

if (save == 1) && (show == 1)
    fig = figure("Name", "Polarization Images", 'Visible','on');
elseif (show == 1) && (save==0)
    fig = figure("Name", "Polarization Images", 'Visible','on');
elseif (save == 1 ) && (show == 0)
    fig = figure("Name", "Polarization Images", 'Visible','off');
end

ax1 = subplot(221);
imshow(demosaic(uint8(proc1), 'rggb'), []);
title("90°")

ax2 = subplot(222);
imshow(demosaic(uint8(proc2), 'rggb'), []);
title("45°")

ax3 = subplot(223);
imshow(demosaic(uint8(proc3), 'rggb'), []);
title("135°")

ax4 = subplot(224);
imshow(demosaic(uint8(proc4), 'rggb'), []);
title("0°")

    if save == 1
        exportgraphics(ax1, fullfile(output_path, [baseFilename '_90.png']), 'Resolution', 400);
        exportgraphics(ax2, fullfile(output_path, [baseFilename '_45.png']), 'Resolution', 400);
        exportgraphics(ax3, fullfile(output_path, [baseFilename '_135.png']), 'Resolution', 400);
        exportgraphics(ax4, fullfile(output_path, [baseFilename '_0.png']), 'Resolution', 400);

        exportgraphics(fig, fullfile(output_path, [baseFilename '_all.png']), 'Resolution', 400);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pol_proc()
    show = 1; %change to 1 if you want to see the RAW image and the separate polarization images
    save = 1; %chage to save or not to save
    input_file = '../images/test/krystal4.raw';
    output_path = './output_images/krystal4/';
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    Z = load_image(show, output_path, input_file, save); %change as needed
    demosaic_polarization_image(Z,output_path, input_file, save, show);
    pause(0.1)
    [DoLP, AoLP] = calculate_polarization(Z);
    visualize_polarization(DoLP, AoLP, output_path, input_file, save);
end

pol_proc;