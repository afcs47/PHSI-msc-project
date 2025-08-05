function [validROI, pos_rect] = select_roi(hsfiltered,R,G,B)
    f = figure('WindowState', 'maximized');
    imshow(hsfiltered(:,:,[R, G, B]),[]) % Displays the image for ROI selection
    title('Select ROI:')
    h_rect = imrect(); % Allows the user to select the ROI
    pos_rect = round(h_rect.getPosition()); % Gets and rounds the ROI position
    close(f) % Closes the figure

    % Check if the ROI is within bounds
    if pos_rect(2) + pos_rect(4) <= size(hsfiltered,1)... %% Check height limit
        && pos_rect(1) + pos_rect(3) <= size(hsfiltered,2)... % Check width limit
        && pos_rect(2) >= 1 && pos_rect(1) >= 1 % Check non-negative index: top-left corner is not above row 1 and top-left corner is not left of column 1
        
        validROI = true; % Exit loop if ROI is valid
    else
        warndlg('Selected ROI is out of bounds. Please select again.', 'Invalid ROI');
        validROI = false;
    end
end
