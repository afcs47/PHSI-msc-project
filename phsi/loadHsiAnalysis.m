function [output1, output2, wavelengths] = loadHsiAnalysis(method, selectedType, datasetName) % process HSI data according to selected method ('Standard', 'Fourier','SPIE', 'SPIE Simplified' or 'Reflectance'), 
% returning output1 - DoLP map (for polarization method) or reflectances (for 'Reflectance') and output2 - AoLP map or []
    
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

        case 'reflectance'
            % Just return the resized reflectances
            reflectancesFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('Spatial_Reflectances_Results_%s*.mat', selectedType)));
            if isempty(reflectancesFile); error('Error: HSI reflectances file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile(['hsi results+figures\' datasetName], reflectancesFile(end).name), 'Spatial_reflectances', 'wavelengths');
            output1 = Spatial_reflectances;
            output2 = [];
        
        otherwise
            error('Error: Unknown HSI method: %s', method);
    end

end
