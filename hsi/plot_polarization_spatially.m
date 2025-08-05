
function plot_polarization_spatially(data, colormap_name, titleText)   
    % Compute common color scale limits
    %Lim = [min(data(:)), max(data(:))]/max(data(:)); 
    if contains(titleText, 'DoLP')
        Lim = [0, 1]; %DoLP
    else; Lim = [-90, 90]; end %AoLP

    fig = figure('Name', titleText);
    %imagesc(data);
    imshow(data, [])
    colorbar;
    colormap(colormap_name);
    title(titleText);
    xlabel('X Position'); 
    ylabel('Y Position');
    axis image; % Ensure square pixels
    clim(Lim); % Set a common color scale

    export_figure(fig, [strrep(titleText, ' ', '_')], 'figures');

end
