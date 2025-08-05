%{
load_image - serves the purpose of loading the image in a usable format for
further processing (uint8). works only for specified dimensions as that is 
the format of the sensor used in our camera. Incorporates substraction of
Dark for calibration

@input_file - is path to file specified for processing
@output_path - path for data exporting
@show - 1=show, 0=don't show
@save - 1=save, 0=don't save
%}

function Z = load_image(input_file, output_path, show, save)
    [~, baseFilename, ~] = fileparts(input_file);
    row = 2048;  col = 2448;

    % load dark for calibration
    dark_path = fopen("dark.raw","r");
    D = fread(dark_path, row * col, 'uint8=>uint8');
    Dark = reshape(D, col, row);
    Dark = Dark';

    % load image with given dimensions
    
    file_in = fopen(input_file, 'r');
    I = fread(file_in, row * col, 'uint8=>uint8'); 
    Z = reshape(I, col, row);
    Z = Z';

    Z = Z - Dark;
        if show == true 
            figure('Visible','on')
            imshow(Z);
            ax = gca;
        else
            figure(Visible="off")
            imshow(Z);
            ax = gca;
            % title("Input RAW image") %uncomment as is necessary
        end
        if save == 1 
            exportgraphics(ax, fullfile(output_path, [baseFilename '_raw.png']), 'Resolution', 400);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
calculate_polarization - is used as a helper function to calculate DoLP and
AoLP, to reduce artifacts in the calculated image we demosaic each set of
polarized pixels separately. Then we calculate the Stokes vector from given
intensities and we use those to calculate DoLP and AoLP for each pixel

@Z - raw image to use for processing
%}

function [DoLP, AoLP] = calculate_polarization(Z)
    % extract the four polarization orientations
    pol_90 = Z(1:2:end, 1:2:end);  % top-left pixels
    pol_45 = Z(1:2:end, 2:2:end);  % top-right pixels
    pol_135 = Z(2:2:end, 1:2:end); % bottom-left pixels
    pol_0 = Z(2:2:end, 2:2:end);   % bottom-right pixels
    
    % demosaic each polarization separately
    color_90 = demosaic(uint8(pol_90), 'rggb');
    color_45 = demosaic(uint8(pol_45), 'rggb');
    color_135 = demosaic(uint8(pol_135), 'rggb');
    color_0 = demosaic(uint8(pol_0), 'rggb');
    
    [height, width, channels] = size(color_0);
    
    % init DoLP and AoLP
    DoLP = zeros(height, width, channels);
    AoLP = zeros(height, width, channels);
    
    % process each color channel
    for c = 1:channels
        % intensities for each polarization angle (need to be type double to
        % work in sqrt function)
        I_0 = double(color_0(:,:,c));
        I_90 = double(color_90(:,:,c));
        I_45 = double(color_45(:,:,c));
        I_135 = double(color_135(:,:,c));
        
        %Stokes vector
        S0 = I_0 + I_90;  % total image intensity
        S1 = I_0 - I_90;  % diff between horizontal and vertical polarization
        S2 = I_45 - I_135; % diff between 45° and 135° polarization
        
        DoLP(:,:,c) = sqrt(S1.^2 + S2.^2) ./ (S0 + eps); %eps to not divide by zero
        
        AoLP(:,:,c) = 0.5 * atan2(S2, S1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
visualize_polarization - serves to output figures of DoLP and AoLP in a
similar format as to what is usually used. 
For DoLP we calculate a mean of
all the color channels and display it with the colormap jet, higher value,
more polarization, is red, lower value, less polarization, is blue
AoLP is first made into x and y coordinates (because of the angle
periodicity) then we calculate the angle from these x,y coordinates back 
again. We use HSV to set Hue according to the angles, remapping the angles 
to [0,1] to use the whole colour spectrum in the angle range. Saturation is 
always set to 1, Value is set to 0 outside the angle range

@DoLP - Degree of Linear Polarization
@AoLP - Angle of Linear Polarization
@input_file - original input file path, for base filename
@output_path - path for data exporting
@angle_min - smallest angle of AoLP to work with
@angle_max - biggest angle of AoLP to work with
@show - 1=show, 0=don't show
@save - 1=save, 0=don't save
%}

function visualize_polarization(DoLP, AoLP, input_file, output_path, angle_min, angle_max, show, show_aolp, save, save_aolp)
    angle_min = angle_min * (pi/180);
    angle_max = angle_max * (pi/180);

    [~, baseFilename, ~] = fileparts(input_file);
    [height, width, channels] = size(DoLP);
    
    %DoLP figure
    if show == 0
        dolp_fig = figure('Name', 'Mean DoLP', 'Visible', 'off');
    else
        dolp_fig = figure('Name', 'Mean DoLP', 'Visible', 'on');
    end
    meanDoLP = mean(DoLP, 3);
    imshow(meanDoLP);
    colormap(jet);
    colorbar;
    
    % AoLP figure
    if show_aolp == 0
        aolp_fig = figure('Name', 'Mean AoLP', 'Visible', 'off');
    else
        aolp_fig = figure('Name', 'Mean AoLP', 'Visible', 'on');
    end
    
    % axes for AoLP so there is space for reference colorwheel
    % aolp_ax = axes('Position', [0.05, 0.05, 0.75, 0.9]);
    aolp_ax = axes;
    
    % init x,y coordinates
    x_total = zeros(height, width);
    y_total = zeros(height, width);
    
    for c = 1:channels
        % convert to x,y because of the periodicity of AoLP
        x_component = cos(2 * AoLP(:,:,c));
        y_component = sin(2 * AoLP(:,:,c));
        
        % sum of all components
        x_total = x_total + x_component;
        y_total = y_total + y_component;
    end
    
    % convert back to angle from x,y
    mean_AoLP = 0.5 * atan2(y_total, x_total);
    mean_AoLP = mod(mean_AoLP, pi);
    
    % mask for angles withing the angle range
    angle_mask = (mean_AoLP >= angle_min) & (mean_AoLP <= angle_max);
    
    % init hsv image
    hsv_img = zeros(height, width, 3);
    
    % remap angles from angle range to [0,1] to use whole color spectrum
    angle_range = angle_max - angle_min;
    if angle_range > 0
        normalized_angles = (mean_AoLP - angle_min) / angle_range;
        hsv_img(:,:,1) = normalized_angles;
    else
        hsv_img(:,:,1) = 0;
    end
    
    hsv_img(:,:,2) = ones(size(mean_AoLP)); % full saturation
    
    % set value to 0 for pixels out of the angle range according to the mask
    hsv_img(:,:,3) = angle_mask;
    
    rgb_img = hsv2rgb(hsv_img);
    imshow(rgb_img, 'Parent', aolp_ax);
    
    % modify colorbar for angle referencing
    polarization_colorbar(aolp_ax, angle_min, angle_max, angle_range);
    
    if save == 1
        exportgraphics(dolp_fig, fullfile(output_path, [baseFilename '_dolp.png']), 'Resolution', 400);
    end
    if save_aolp == 1
        exportgraphics(aolp_fig, fullfile(output_path, [baseFilename '_aolp.png']), 'Resolution', 400);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
polarization_colorbar - customizes a standard matlab colorbar to have the
same colormap as the AoLP image and shows ten linearly spaced tick marks to
make reading out the angle easier

@parent_ax - AoLP image axes, where to put colorbar
@angle_min - smallest angle of AoLP to work with
@angle_max - biggest angle of AoLP to work with
@angle_range - difference between max and min angles
%}

function polarization_colorbar(parent_ax, angle_min, angle_max, angle_range)    
    cb = colorbar(parent_ax);

    n_points = 256;   
    angles = linspace(angle_min, angle_max, n_points)';
    
    % create HSV color map
    hsv_map = zeros(n_points, 3);
    
    % Map angles to hue values [0,1]
    normalized_angles = (angles - angle_min) / angle_range;
    hsv_map(:,1) = normalized_angles;
    hsv_map(:,2) = 1;
    hsv_map(:,3) = 1;
    
    rgb_map = hsv2rgb(hsv_map);
    
    % colormap the parent axes
    colormap(parent_ax, rgb_map);
    
    % limits of the parent axes
    clim(parent_ax, [angle_min, angle_max]);
    
    % linearly spaced tick marks
    num_ticks = 10;
    tick_angles = linspace(angle_min, angle_max, num_ticks);
    tick_degrees = round(tick_angles * (180/pi),2,"decimals");
    
    % print the tick marks and labels
    cb.Ticks = tick_angles;
    cb.TickLabels = arrayfun(@(x) sprintf('%g°', x), tick_degrees, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
demosaic_polarization_image - function to display only the pixels with
given polarizer orientations

@Z - raw image to use for processing
@input_file - original input file path, for base filename
@output_path - path for data exporting
@show - 1=show, 0=don't show
@save - 1=save, 0=don't save
%}

function demosaic_polarization_image(Z, input_file, output_path, show, save)
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

    %hide if you dont need full figure
    exportgraphics(fig, fullfile(output_path, [baseFilename '_all.png']), 'Resolution', 400);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
pol_proc - main function used for calling all the others, user specifies if
they want to save/ show the figures or both, input file to process, output
path where the potentionally saved images will be found and angle range
that user wants displayed
%}
function main()
    show = 1; %change to 1 if you want to see the RAW image and the separate polarization images
    save = 1; %change to save or not to save only the polarization images
    show_aolp = 1; % show only polarimetry (AoLP) image
    save_aolp = 0; % save only polarimetry (AoLP) image

    %change path to file accordingly
    input_file = '../images/roztoky/laser.raw';
    output_path = '../output_images/roztoky/laser/'; % for all images
    % output_path = '../output_images/roztoky/aolp_calcite/'; % only aolp
    
    %min and max angle (degrees)
    angle_min = 0; 
    angle_max = 180; % max 180

    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    Z = load_image(input_file, output_path, show, save);
    demosaic_polarization_image(Z, input_file, output_path, show, save);
    pause(0.2)
    [DoLP, AoLP] = calculate_polarization(Z);
    visualize_polarization(DoLP, AoLP, input_file, output_path, angle_min, angle_max, show, show_aolp, save, save_aolp);
end

main;