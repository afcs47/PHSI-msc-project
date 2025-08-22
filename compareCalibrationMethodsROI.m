function [tform_best, tform_fit, tform_est, report] = compareCalibrationMethodsROI(ptsHSI, ptsPOL, imagePointsHSI, imagePointsPOL, hsi_rgb, pol_rgb, check_calibration) % Compare fitgeotform2d (4 manual corners) vs estgeotform2d (matched ROI
% points) from the input of ptsHSI, ptsPOL (4x2 corners selected from detected points), imagePoints1 (HSI), imagePoints2 (POL) (all detected
% points from detectCheckerboardPoints) and return the better transform for POL->HSI alignment, as well as a struct variable with errors and counts

    %% Corner-only transform - fitgeotform2d
    tform_fit = fitgeotform2d(ptsPOL, ptsHSI, "projective"); % Fit a projective transform
    [xf, yf] = transformPointsForward(tform_fit, ptsPOL(:,1), ptsPOL(:,2)); % Applies the transform to POL corners to see where they map in the HSI image
    err_fit = sqrt(mean(sum(([xf yf] - ptsHSI).^2, 2))); % RMSE between warped POL corners and HSI corners

    %% ROI points 
    % Finds which detected points are contained within the quadrilateral defined by the 4 manually selected corners

    % Organize points into grid-ish order: sort by x (cols), then y (rows) ie starting in the most TL point and ending in the most BR point - get the points indexes and for plotting later
    [~, idxHSI] = sortrows(imagePointsHSI, [1,2]); % column-major ordering (top-to-bottom, then move right)
    sortedHSI = imagePointsHSI(idxHSI,:);
    [~, idxPOL] = sortrows(imagePointsPOL, [1,2]);
    sortedPOL = imagePointsPOL(idxPOL,:);

    % Identify the indexes of the selected corners
    idxPtsHSI = findIndexOfSelectedPoints(ptsHSI, imagePointsHSI); % row-major ordering (left-to-right, then move down)
    idxPtsPOL = findIndexOfSelectedPoints(ptsPOL, imagePointsPOL);

    % Determine which checkerboard is smaller
    nHSI = size(imagePointsHSI,1);
    nPOL = size(imagePointsPOL,1);
    
    % Find bounding box of smaller checkerboard and cropp points from larger one inside this ROI
    if nHSI <= nPOL
      [roiHSI, roiPOL] = selectCommonPoints(sortedHSI, sortedPOL, idxPtsHSI, idxPtsPOL);
    else
      [roiPOL, roiHSI] = selectCommonPoints(sortedPOL, sortedHSI, idxPtsPOL, idxPtsHSI);
    end

    % Enforce it only keeps the minimum/common nr of points
    nKeep = min(size(roiHSI,1), size(roiPOL,1));
    roiHSI = roiHSI(1:nKeep,:);
    roiPOL = roiPOL(1:nKeep,:);
    
    fprintf('\nSelected %d common points (min of HSI=%d, POL=%d).\n', nKeep, size(roiHSI,1), size(roiPOL,1));

    %% More-points transform - estgeotform2d
    tform_est = estgeotform2d(roiPOL, roiHSI, "projective"); % Estimate a projective transform
    [xe, ye] = transformPointsForward(tform_est, roiPOL(:,1), roiPOL(:,2)); % Applies the transform to POL corners to see where they map in the HSI image
    err_est = sqrt(mean(sum(([xe ye] - roiHSI).^2, 2))); % RMSE between warped POL points and HSI points

    %% Choose best result by comparing RMSE of both methods
    tform_best = tform_fit;
    if err_est < err_fit
        tform_best = tform_est;
    end

    if check_calibration
        %% Visualizations
        % Overlay plot
        figure('Name','Warped POL Points (red = fit, blue = est)');
        imshow(hsi_rgb); axis image; hold on;
        plot(roiHSI(:,1), roiHSI(:,2), 'g.', 'MarkerSize', 12); %green squares - detected HSI points (ROI + 4 corners)
        plot(ptsHSI(:,1), ptsHSI(:,2), 'go', 'MarkerSize', 10, 'LineWidth', 2); %green circler - 4 corners
        plot(xf, yf, 'rx', 'MarkerSize', 12, 'LineWidth', 2); %red x - POL points warped using 4 corner fit
        plot(xe, ye, 'b+', 'MarkerSize', 8, 'LineWidth', 1.5); %blue + - POL points warped using ROI (more points) fit
        legend('HSI ROI points','HSI 4 corners','fitgeotform2d warped','estgeotform2d warped');
        title(sprintf('RMSE fit=%.2f px, est=%.2f px (best: %s)', err_fit, err_est, ternary(err_est<err_fit,'est','fit')));
        hold off;
    
        % Side-by-side point match plot
        figure('Name','Matched ROI points side-by-side');
        subplot(1,2,1);
        imshow(pol_rgb); axis image; hold on;
        plot(roiPOL(:,1), roiPOL(:,2), 'r.', 'MarkerSize', 10);
        plot(sortedPOL(:,1), sortedPOL(:,2), 'ro', 'MarkerSize', 6);
        title('POL ROI points (red = matched)');
        subplot(1,2,2);
        imshow(hsi_rgb); axis image; hold on;
        plot(roiHSI(:,1), roiHSI(:,2), 'g.', 'MarkerSize', 10);
        plot(sortedHSI(:,1), sortedHSI(:,2), 'go', 'MarkerSize', 6);
        title('HSI ROI points (green = matched)');
    
        %Compare alignment methods side-by-side
        alignedPOL_fit = imwarp(pol_rgb, tform_fit, 'OutputView', imref2d(size(hsi_rgb)));
        alignedPOL_est = [];
        if ~isempty(tform_est)
            alignedPOL_est = imwarp(pol_rgb, tform_est, 'OutputView', imref2d(size(hsi_rgb)));
        end
        
        figure('Name','Compare alignment methods side-by-side');
        title('Compare alignment methods side-by-side');
        subplot(1,2,1); imshowpair(hsi_rgb, alignedPOL_fit, 'falsecolor'); title('fitgeotform2d result'); axis image;
        if ~isempty(alignedPOL_est)
            subplot(1,2,2); imshowpair(hsi_rgb, alignedPOL_est, 'falsecolor'); title('estgeotform2d result'); axis image;
        end
    end

    %% Report
    report.err_fit = err_fit;
    report.err_est = err_est;
    report.n_roi_hsi = size(roiHSI,1);
    report.n_roi_pol = size(roiPOL,1);
    
    if check_calibration; disp(report); end
end

% Small inline ternary operation (text in figure's title)
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

% Normalize function
function ptsNorm = normalizePoints(pts)
    % Scale x,y independently to [0,1]
    minVals = min(pts,[],1);
    maxVals = max(pts,[],1);
    ptsNorm = (pts - minVals) ./ (maxVals - minVals);
end

function [roiRef, roiCmp] = selectCommonPoints(imagePointsRef, imagePointsCmp, idxPtsRef, idxPtsCmp) %ref =  smaller checkerboard points (reference), cmp = larger checkerboard points (compared)
    % Normalize coordinates to [0,1] for easier comparison
    refNorm = normalizePoints(imagePointsRef);
    cmpNorm = normalizePoints(imagePointsCmp);

    % Build masks of points inside ROI
    maskCmp = getROIMask(cmpNorm, refNorm, idxPtsCmp);
    maskRef = getROIMask(refNorm, refNorm, idxPtsRef);

    % Keep only ROI points
    roiCmp = imagePointsCmp(maskCmp,:);
    roiRef = imagePointsRef(maskRef,:);
end

function idxPts = findIndexOfSelectedPoints(pts, imagePoints)
    % Find indices of the selected HSI/POL points in the respective imagePoints...
    idxPts = zeros(size(pts,1),1);
    for k = 1:size(pts,1)
        idxPts(k) = find(ismember(imagePoints, pts(k,:), 'rows'));
    end
end

function mask = getROIMask(norm, refNorm, idxPts)
    %Extract corners from normalized points
    corners = [norm(idxPts(1),:);  % top left
           norm(idxPts(2),:);  % top right
           norm(idxPts(3),:);  % bottom right
           norm(idxPts(4),:)]; % bottom left

    % Define grid size from refNorm
    nColsRef = numel(unique(refNorm(:,1))); % points along x
    nRowsRef = numel(unique(refNorm(:,2))); % points along y

    % Normalized parametric grid - [0,1]x[0,1]
    [u,v] = meshgrid(linspace(0,1,nColsRef), linspace(0,1,nRowsRef));

    % Bilinear interpolation of corners (inside quadrilateral) - expected grid in input space
    gridX = (1-u).*((1-v)*corners(1,1) + v*corners(4,1)) + ...
             u.*((1-v)*corners(2,1) + v*corners(3,1));
    gridY = (1-u).*((1-v)*corners(1,2) + v*corners(4,2)) + ...
             u.*((1-v)*corners(2,2) + v*corners(3,2));
    
    gridPts = [gridX(:), gridY(:)];

    % Find nearest-neighbor of cmpNorm points for each grid location
    idxMatched = knnsearch(norm, gridPts);

    % Mask selecting exactly those points
    mask = false(size(norm,1),1);
    mask(idxMatched) = true;
end