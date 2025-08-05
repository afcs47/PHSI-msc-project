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
