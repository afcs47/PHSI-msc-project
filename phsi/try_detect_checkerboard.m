function [imagePoints, boardSize] = try_detect_checkerboard(img, camType) % Attempt Checkerboard Detection through multiple enhancement techniques, using Computer Vision Toolbox's predefined function
    
    %size(img)

    % Try original image
    [imagePoints, boardSize] = detectCheckerboardPoints(img);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (original/ no enhancement needed).\n", camType);
        return;
    end

    % Try grayscale
    if size(img,3) == 3
        gray = rgb2gray(img); %"MAP must be a m x 3 array"
    else
        gray = im2gray(img); %"Use im2gray for RGB and grayscale images"
    end
    %size(gray)
    [imagePoints, boardSize] = detectCheckerboardPoints(gray);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (grayscale).\n", camType);
        return;
    end

    % Try contrast enhancement
    enhanced = imadjust(gray);
    %size(enhanced)
    [imagePoints, boardSize] = detectCheckerboardPoints(enhanced);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (enhanced contrast).\n", camType);
        return;
    end

    % Try edge sharpening
    sharp = imsharpen(enhanced);
    %size(sharp)
    [imagePoints, boardSize] = detectCheckerboardPoints(sharp);
    if ~isempty(imagePoints)
        fprintf("✅ Checkerboard detected in %s image (sharpened).\n", camType);
    else
        fprintf("⚠️ Checkerboard NOT detected in %s image after all attempts.\n", camType);
    end
end