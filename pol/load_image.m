%{
load_image - serves the purpose of loading the image in a usable format for
further processing (uint8). works only for specified dimensions as that is 
the format of the sensor used in our camera. Incorporates substraction of
Dark for calibration

@input_file - is path to file specified for processing
@output_path - path for data exporting
@show - 1=show, 0=don't show
@save - 1=save, 0=don't save
%}

function Z = load_image(input_file, output_path, show, save, dark_path) % Adapted from "polarimetry.m" to only load and subtract dark image if dark_path is non-empty; handle error and file closing; simplify image display; added figure name
    [~, baseFilename, ~] = fileparts(input_file);
    row = 2048;  col = 2448;

    % load image with given dimensions   
    file_in = fopen(input_file, 'r');
    I = fread(file_in, row * col, 'uint8=>uint8'); 
    fclose(file_in);
    
    Z = reshape(I, col, row); % Reshape the linear array into a 2D matrix ([col x row])
    Z = Z'; % Transpose the image to correct the orientation to [row x col]

    % Apply dark calibration only if path is provided
    if nargin >= 5 && ~isempty(dark_path)
        dark_file = fopen(dark_path, 'r'); % load dark for calibration
        if dark_file == -1 % handle error
            error('Could not open dark file at path: %s', dark_path);
        end
        D = fread(dark_file, row * col, 'uint8=>uint8');
        fclose(dark_file);

        Dark = reshape(D, col, row);
        Dark = Dark';

        %imshow(Z);
        %title("Dark RAW image")

        Z = Z - Dark;
    end

    
    % Display image
    if show
        figure("Name","Input RAW image", 'Visible', 'on');
    else
        figure("Name","Input RAW image", 'Visible', 'off');
    end
    imshow(Z);
    ax = gca;
    title("Input RAW image")

    % Save image if required
    if save
        exportgraphics(ax, fullfile(output_path, [baseFilename '_raw.png']), 'Resolution', 400);
    end
end
