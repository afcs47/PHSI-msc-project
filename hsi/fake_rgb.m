function [R, G, B] = fake_rgb(showRGB, HS_calibrated, datasetname)
    binning = 8; % Set original binning factor, used for downsampling
    
    % Calculate the approximate indices for the RGB channels
    R = round(333/binning); % channel number (not wavelenght) / binning factor;
    G = round(205/binning); % channel number (not wavelenght) / binning factor;
    B = round(100/binning); % channel number (not wavelenght) / binning factor;
    
    if showRGB
        % Display the calibrated data as a pseudo-RGB image
        figure();
        imshow(HS_calibrated(:,:,[R, G, B]),[]);
        title(['Pseudo-RGB for' datasetname]);
    end 
end

    