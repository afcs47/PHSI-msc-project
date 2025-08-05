
function mean_reflectances = resize_reflectances(mean_reflectances)

    % Get all field names (e.g., angle_0, angle_45, ...)
    fields = fieldnames(mean_reflectances);
    
    % Determine the smallest common size
    minRows = inf;
    minCols = inf;
    
    for i = 1:numel(fields)
        [rows, cols] = size(mean_reflectances.(fields{i}));
        minRows = min(minRows, rows);
        minCols = min(minCols, cols);
    end
    
    % Crop each image to the smallest size
    for i = 1:numel(fields)
        data = mean_reflectances.(fields{i});
        mean_reflectances.(fields{i}) = data(1:minRows, 1:minCols);  % Crop to min size
    end

end