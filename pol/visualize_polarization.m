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

function visualize_polarization(DoLP, AoLP, input_file, output_path, angle_min, angle_max, show, show_aolp, save, save_aolp, method) % added 'method' input
    angle_min = angle_min * (pi/180);
    angle_max = angle_max * (pi/180);

    [~, baseFilename, ~] = fileparts(input_file);
    [height, width, channels] = size(DoLP);
    
    %DoLP figure
    if show == 0
        dolp_fig = figure('Name', ['Mean DoLP - ' method], 'Visible', 'off');
    else
        dolp_fig = figure('Name', ['Mean DoLP - ' method], 'Visible', 'on');
    end
    meanDoLP = mean(DoLP, 3);
    imshow(meanDoLP);
    colormap(jet);
    colorbar;
    title(['Mean DoLP - ' method]);
    
    % AoLP figure
    if show_aolp == 0
        aolp_fig = figure('Name', ['Mean AoLP - ' method], 'Visible', 'off');
    else
        aolp_fig = figure('Name', ['Mean AoLP - ' method], 'Visible', 'on');
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
    title(['Mean AoLP - ' method]);

    % modify colorbar for angle referencing
    polarization_colorbar(aolp_ax, angle_min, angle_max, angle_range);
    
    if save == 1
        exportgraphics(dolp_fig, fullfile(output_path, [baseFilename regexprep(method, ' ', '_') '_dolp.png']), 'Resolution', 400);
    end
    if save_aolp == 1
        exportgraphics(aolp_fig, fullfile(output_path, [baseFilename regexprep(method, ' ', '_') '_aolp.png']), 'Resolution', 400);
    end

end