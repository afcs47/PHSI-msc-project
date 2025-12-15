function offsetY = estimate_offsetY(refImage, warpedImage) % Interactive Y-offset fine tuning user estimation after homography

    disp('Use ↑ ↓ arrows to adjust the vertical alignment.');
    disp('Press ENTER to confirm, ESC to cancel (returns 0).');

    %refImage = mat2gray(refImage);
    %warpedImage = mat2gray(warpedImage);

     % Create figure and show initial overlay
    fig = figure('Name','Estimate offsetY','NumberTitle','off','KeyPressFcn',@keyHandler);
    % Use false color overlay with transparency
    ax = axes('Parent', fig);
    imshow(refImage, [0 1], 'Parent', ax); colormap(ax, 'jet'); axis image; hold(ax, 'on');
    h = imshow(warpedImage, [0 1], 'Parent', ax);
    set(h, 'AlphaData', 0.65); % Semi-transparent overlay
    title(ax, 'Use ↑ ↓ to align. Press ENTER to confirm.');
    axis(ax, 'image');
    

    offsetY = 0;
    step = 1; % pixel per press

    %Nested callback
    function keyHandler(~, event)
        switch event.Key
            case 'uparrow'
                offsetY = offsetY - step;
            case 'downarrow'
                offsetY = offsetY + step;
            case 'return'
                uiresume(fig);
            case 'escape'
                offsetY = 0;
                uiresume(fig);
        end
        shifted = imtranslate(warpedImage, [0 offsetY]);
        set(h, 'CData', shifted);
        title(ax, sprintf('offsetY = %d px | ↑↓ adjust, ENTER confirm, ESC cancel', offsetY));
        drawnow;
    end

    % Wait for user
    uiwait(fig);
    if isvalid(fig), close(fig); end

    %disp(['Final estimated offsetY = ', num2str(offsetY)]);
end
