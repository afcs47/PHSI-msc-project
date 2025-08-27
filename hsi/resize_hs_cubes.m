function cube = resize_hs_cubes(hs_cubes) %similar to resize_reflectances.m but for 3d datasets

    fields = fieldnames(hs_cubes);

    % Determine the smallest common spatial size
    minRows = inf; minCols = inf;
    for i = 1:numel(fields)
        [rows, cols, ~] = size(hs_cubes.(fields{i}));
        minRows = min(minRows, rows);
        minCols = min(minCols, cols);
    end

    % Initialize running sum
    first = hs_cubes.(fields{1});
    cubeSum = zeros(minRows, minCols, size(first,3), 'like', first);

    % Crop each cube to the smallest size
    for i = 1:numel(fields)
        data = hs_cubes.(fields{i});
        %hs_cubes.(fields{i}) = data(1:minRows, 1:minCols, :);
        data = data(1:minRows, 1:minCols, :); % crop
        %hs_cubes.(fields{i}) = data; % assign cropped
        cubeSum = cubeSum + data(1:minRows, 1:minCols, :); % Sum directly over cropped cubes
        clear data; % free immediately to avoid 'Out of memory.' error
    end
    cube = cubeSum ./ numel(fields);
end
