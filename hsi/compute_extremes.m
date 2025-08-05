
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