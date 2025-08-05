
%{
pol_proc - main function used for calling all the others, user specifies if
they want to save/ show the figures or both, input file to process, output
path where the potentionally saved images will be found and angle range
that user wants displayed
%}
function pol_proc() %also called main in "polarimetry.m"
    show = 1; %change to 1 if you want to see the RAW image and the separate polarization images
    save = 1; %change to save or not to save only the polarization images
    show_aolp = 1; % show only polarimetry (AoLP) image
    save_aolp = 0; % save only polarimetry (AoLP) image

    %change path to file accordingly
    input_file = '../images/roztoky/laser.raw';
    output_path = '../output_images/roztoky/laser/'; % for all images
    % output_path = '../output_images/roztoky/aolp_calcite/'; % only aolp
    
    %min and max angle (degrees)
    angle_min = 0; 
    angle_max = 180; % max 180

    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    Z = load_image(input_file, output_path, show, save);
    demosaic_polarization_image(Z, input_file, output_path, show, save);
    pause(0.2)
    [DoLP, AoLP] = calculate_polarization(Z);
    visualize_polarization(DoLP, AoLP, input_file, output_path, angle_min, angle_max, show, show_aolp, save, save_aolp);
end
