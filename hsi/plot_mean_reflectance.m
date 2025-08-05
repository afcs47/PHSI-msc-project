function plot_mean_reflectance(datasetname, wavelengths, mean_ref, std_ref, choice_std_dev, choice_fill, colors, k)
    plot(wavelengths, mean_ref,'Color', colors(k,:), 'DisplayName', datasetname); % Plots the mean reflectance for the given dataset      
    if strcmp(choice_std_dev, 'Yes')
        if strcmp(choice_fill, 'No')
            plot(wavelengths, mean_ref + std_ref, '-', 'Color', colors(k,:), 'DisplayName', strcat(datasetname, ' Max')); % Plots the upper bound (mean + std)
            plot(wavelengths, mean_ref - std_ref, '-', 'Color', colors(k,:), 'DisplayName', strcat(datasetname, ' Min')); % Plots the lower bound (mean - std)
        elseif strcmp(choice_fill, 'Yes')
            % Plot mean reflectance with shaded standard deviation
            x_vals = [wavelengths, fliplr(wavelengths)]; % X values for fill (wavelengths mirrored)
            y_vals = [mean_ref - std_ref, fliplr(mean_ref + std_ref)]; % Y values for fill
            
            % Plot shaded region for standard deviation
            fill(x_vals, y_vals, colors(k,:), 'FaceAlpha', 0.2, 'EdgeColor', colors(k,:), 'DisplayName', strcat(datasetname, ' Std Dev')); 
        end
    end 
end
