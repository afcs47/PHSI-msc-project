function [imagePoints, boardSize] = tryDetectCheckerboard(img, camType) % Attempt Checkerboard Detection through multiple enhancement techniques, using Computer Vision Toolbox's predefined function
    % Try original image
    [imagePoints, boardSize] = detectCheckerboardPoints(img);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (original/ no enhancement needed).\n", camType);
        return;
    end

    % Try grayscale
    gray = rgb2gray(img);
    [imagePoints, boardSize] = detectCheckerboardPoints(gray);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (grayscale).\n", camType);
        return;
    end

    % Try contrast enhancement
    enhanced = imadjust(gray);
    [imagePoints, boardSize] = detectCheckerboardPoints(enhanced);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (enhanced contrast).\n", camType);
        return;
    end

    % Try edge sharpening
    sharp = imsharpen(enhanced);
    [imagePoints, boardSize] = detectCheckerboardPoints(sharp);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (sharpened).\n", camType);
    else
        fprintf("⚠️ Checkerboard NOT detected in %s image after all attempts.\n", camType);
    end
end