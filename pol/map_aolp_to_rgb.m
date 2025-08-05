function rgb_img = map_aolp_to_rgb(mean_AoLP, angle_min, angle_max) %Based on visualize_polarization function from polarimetry.m
    % Ensure angles are in radians
    angle_min = deg2rad(angle_min);
    angle_max = deg2rad(angle_max);

    % Mask for valid angles
    angle_mask = (mean_AoLP >= angle_min) & (mean_AoLP <= angle_max);
    angle_range = angle_max - angle_min;

    % Normalize angle to [0,1] hue
    normalized_angles = (mean_AoLP - angle_min) / angle_range;

    % HSV image init
    hsv_img = zeros([size(mean_AoLP), 3]);
    hsv_img(:,:,1) = normalized_angles;
    hsv_img(:,:,2) = 1;
    hsv_img(:,:,3) = double(angle_mask);  % Brightness mask

    rgb_img = hsv2rgb(hsv_img);
end