function [output1, output2] = load_hsi_analysis(method, selectedType, datasetName, reflectanceMode) % process HSI data according to selected method ('Standard', 'FFT', 'LSQ' or 'Reflectance'), 
% returning output1 - DoLP map (for polarization method) or reflectances (for 'Reflectance') and output2 - AoLP map or wavelengths
    
    % Handle optional input
    if nargin < 4
        reflectanceMode = 'Spatial_reflectances'; % default
    end
    
    methodLower = lower(strip(method)); % Normalize the method name (lowercase, removed of all whitespace siding the str)

    %% Branch by method
    switch methodLower
        case 'standard'
            hsiDataFile = dir(fullfile(['hsi results+figures/' datasetName], sprintf('Standard_Stokes_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (Standard method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile('hsi results+figures', datasetName, hsiDataFile(end).name), 'DoLP', 'AoLP');
            output1 = DoLP;
            output2 = AoLP;

        case 'fft'
            hsiDataFile = dir(fullfile(['hsi results+figures/' datasetName], sprintf('Fourier_FFT_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (Fourier FFT method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile('hsi results+figures', datasetName, hsiDataFile(end).name), 'DoLP', 'AoLP');
            output1 = DoLP;
            output2 = AoLP;

        case 'lsq'
            hsiDataFile = dir(fullfile(['hsi results+figures/' datasetName], sprintf('Fourier_LSQ_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile); error('Error: (Fourier LSQ method) HSI data file for %s not found! Please run the respective script first or check directories.', datasetName); end
            load(fullfile('hsi results+figures', datasetName, hsiDataFile(end).name), 'DoLP', 'AoLP');
            output1 = DoLP;
            output2 = AoLP;
    

        case 'reflectance' % Choose between spatial or wavelength reflectances
%             reflectancesFile = dir(fullfile(['hsi results+figures/' datasetName], sprintf('Reflectances_%s*.mat', selectedType)));
%             if isempty(reflectancesFile); error('Error: HSI reflectances file for %s not found! Please run the respective script first or check directories.', datasetName); end
%             load(fullfile('hsi results+figures', datasetName, reflectancesFile(end).name), 'Spectral_reflectances', 'Spatial_reflectances', 'wavelengths');
% 
%             switch lower(reflectanceMode)
%                 case 'spatial_reflectances' % Returns reflectance maps over space
%                     output1 = Spatial_reflectances;
%                     output2 = [];
%                 case 'spectral_reflectances' % Returns reflectances averaged across wavelength dimension
%                     output1 = Spectral_reflectances;
%                     output2 = wavelengths;
%                 otherwise
%                     error('Unknown reflectanceMode: %s. Use ''Spatial_reflectances'' or ''Spectral_reflectances''.', reflectanceMode);
%             end


            % Search for BOTH naming conventions: Reflectances_<datasetName> or woWG_Reflectance_Results_<datasetName>
            baseFolder = fullfile('hsi results+figures', datasetName);
        
            patternWG = sprintf('Reflectances_%s*.mat', selectedType);
            patternNoWG = sprintf('woWG_Reflectance_Results_%s*.mat', selectedType);
        
            reflectancesFile = dir(fullfile(baseFolder, patternWG));
    
            % If WG-style naming not found, look for woWG-style
            if isempty(reflectancesFile), reflectancesFile = dir(fullfile(baseFolder, patternNoWG)); end
            if isempty(reflectancesFile), error('Error: HSI reflectances file for %s not found!\nLooked for patterns:\n  %s\n  %s', datasetName, patternWG, patternNoWG); end
    
            fileToLoad = fullfile(baseFolder, reflectancesFile(end).name);
            load(fileToLoad, 'Spectral_reflectances', 'Spatial_reflectances', 'wavelengths');
    
            switch lower(reflectanceMode)
                case 'spatial_reflectances'
                    output1 = Spatial_reflectances;
                    output2 = [];
        
                case 'spectral_reflectances'
                    output1 = Spectral_reflectances;
                    output2 = wavelengths;
        
                otherwise
                    error('Unknown reflectanceMode: %s. Use ''Spatial_reflectances'' or ''Spectral_reflectances''.', reflectanceMode);
            end


        %% in case the choosen method isn't recognized
        otherwise
            error('Error: Unknown HSI method: %s', method);
    end

end
