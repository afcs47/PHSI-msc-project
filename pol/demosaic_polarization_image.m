%{
demosaic_polarization_image - function to display only the pixels with
given polarizer orientations

@Z - raw image to use for processing
@input_file - original input file path, for base filename
@output_path - path for data exporting
@show - 1=show, 0=don't show
@save - 1=save, 0=don't save
%}

function [proc1, proc2, proc3, proc4] = demosaic_polarization_image(Z, input_file, output_path, show, save) %  added output variables
[~, baseFilename, ~] = fileparts(input_file);

% Extract polarization images directly
proc1 = Z(1:2:end, 1:2:end);     % 90° (top-left pixels)
proc2 = Z(1:2:end, 2:2:end);     % 45° (top-right pixels)
proc3 = Z(2:2:end, 1:2:end);     % 135° (bottom-left pixels)
proc4 = Z(2:2:end, 2:2:end);     % 0° (bottom-right pixels)

if (save == 1) && (show == 1)
    fig = figure("Name", "Polarization Images", 'Visible','on');
elseif (show == 1) && (save==0)
    fig = figure("Name", "Polarization Images", 'Visible','on');
elseif (save == 1 ) && (show == 0)
    fig = figure("Name", "Polarization Images", 'Visible','off');
end

ax1 = subplot(221);
imshow(demosaic(uint8(proc1), 'rggb'), []);
title("90°")

ax2 = subplot(222);
imshow(demosaic(uint8(proc2), 'rggb'), []);
title("45°")

ax3 = subplot(223);
imshow(demosaic(uint8(proc3), 'rggb'), []);
title("135°")

ax4 = subplot(224);
imshow(demosaic(uint8(proc4), 'rggb'), []);
title("0°")

if save == 1
    exportgraphics(ax1, fullfile(output_path, [baseFilename '_90.png'] ), 'Resolution', 400);
    exportgraphics(ax2, fullfile(output_path, [baseFilename '_45.png'] ), 'Resolution', 400);
    exportgraphics(ax3, fullfile(output_path, [baseFilename '_135.png'] ), 'Resolution', 400);
    exportgraphics(ax4, fullfile(output_path, [baseFilename '_0.png'] ), 'Resolution', 400);

    %hide if you dont need full figure
    exportgraphics(fig, fullfile(output_path, [baseFilename '_all.png'] ), 'Resolution', 400);
end
end
