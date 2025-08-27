function [output1, output2] = loadHsiAnalysis(method, selectedType, datasetName, reflectanceMode) % process HSI data according to selected method ('Standard', 'Fourier','SPIE', 'SPIE Simplified' or 'Reflectance'), 
% returning output1 - DoLP map (for polarization method) or reflectances (for 'Reflectance') and output2 - AoLP map or wavelengths
    
    % Handle optional input
    if nargin < 4
        reflectanceMode = 'Spatial_reflectances'; % default
    end

    %% Branch by method
    switch lower(method)
        case 'standard'
            hsiDataFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('Standard_Stokes_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (Standard method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile(['hsi results+figures\' datasetName], hsiDataFile(end).name), 'Standard_DoLP_map', 'Standard_AoLP_map', 'Stokes_S0', 'Stokes_S1', 'Stokes_S2', 'wavelengths');
            output1 = Standard_DoLP_map;
            output2 = Standard_AoLP_map;

        case 'fourier'
            hsiDataFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('Fourier_Stokes_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (Fourier method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile('hsi results+figures', hsiDataFile(end).name), 'Fourier_DoLP_map', 'Fourier_AoLP_map', 'Stokes_S0', 'Stokes_S1', 'Stokes_S2', 'wavelengths');
            output1 = Fourier_DoLP_map;
            output2 = Fourier_AoLP_map;

        case 'spie'
            hsiDataFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('SPIE_Fourier_Stokes_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (SPIE method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile('hsi results+figures', hsiDataFile(end).name), 'SPIE_Fourier_DoLP_map', 'SPIE_Fourier_AoLP_map', 'Stokes_S0', 'Stokes_S1', 'Stokes_S2', 'wavelengths');
            output1 = SPIE_Fourier_DoLP_map;
            output2 = SPIE_Fourier_AoLP_map;

        case 'spie simplified'
            hsiDataFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('SpieSimple_Fourier_Stokes_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (SPIE Simplified method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile('hsi results+figures', hsiDataFile(end).name), 'SpieSimple_Fourier_DoLP_map', 'SpieSimple_Fourier_AoLP_map', 'Stokes_S0', 'Stokes_S1', 'Stokes_S2', 'wavelengths');
            output1 = SpieSimple_Fourier_DoLP_map;
            output2 = SpieSimple_Fourier_AoLP_map;

        case 'reflectance' % Choose between spatial or wavelength reflectances
            reflectancesFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('All_Angles_Reflectances_Results_%s*.mat', selectedType)));
            if isempty(reflectancesFile); error('Error: HSI reflectances file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile(['hsi results+figures\' datasetName], reflectancesFile(end).name), 'Wavelength_reflectances', 'Spatial_reflectances', 'wavelengths');

            switch lower(reflectanceMode)
                case 'spatial_reflectances' % Returns reflectance maps over space
                    output1 = Spatial_reflectances;
                    output2 = [];
                case 'wavelength_reflectances' % Returns reflectances averaged across wavelength dimension
                    output1 = Wavelength_reflectances;
                    output2 = wavelengths;
                otherwise
                    error('Unknown reflectanceMode: %s. Use ''Spatial_reflectances'' or ''Wavelength_reflectances''.', reflectanceMode);
            end
        
        otherwise
            error('Error: Unknown HSI method: %s', method);
    end

end
