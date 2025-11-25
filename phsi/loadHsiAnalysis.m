function [output1, output2] = loadHsiAnalysis(method, selectedType, datasetName, reflectanceMode) % process HSI data according to selected method ('Standard', 'FFT', 'Full FFT', 'LSQ', 'Full LSQ', 'Fourier','SPIE', 'SPIE Simplified' or 'Reflectance'), 
% returning output1 - DoLP map (for polarization method) or reflectances (for 'Reflectance') and output2 - AoLP map or wavelengths
    
    % Handle optional input
    if nargin < 4
        reflectanceMode = 'Spatial_reflectances'; % default
    end
    
    methodLower = lower(strip(method)); % Normalize the method name (lowercase, removed of all whitespace siding the str)

    %% Branch by method
    switch methodLower
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

        %% other methods, saved in the same struct variable and same file
        case {'fft', 'full fft', 'lsq', 'full lsq'}
            % unified structure-based loading
            hsiDataFile = dir(fullfile(['hsi results+figures\' datasetName], sprintf('Fourier_Methods_Results_%s*.mat', selectedType)));
            if isempty(hsiDataFile)
                error('Error: (Fourier methods) HSI file for %s not found!', datasetName);
            end
            %load(fullfile(['hsi results+figures\' datasetName], hsiDataFile(end).name), 'results'); % error: 'results' requires Fixed-Point Designer

            % Load into a temporary struct to avoid name conflicts
            temp = load(fullfile(['hsi results+figures\' datasetName], hsiDataFile(end).name));
            
%             % Extract the actual results structure
%             if isfield(temp, 'results')
%                 dataStruct = temp.results;
%             else
%                 error('The loaded file does not contain a variable named "results".');
%             end
% 
%             varNames = fieldnames(temp);
%             % Attempt to identify the structure variable
%             dataStruct = [];
%             for i = 1:numel(varNames)
%                 v = temp.(varNames{i});
%                 if isstruct(v) && any(isfield(v, {'FFT', 'FullFFT', 'LSQ', 'FullLSQ'}))
%                     dataStruct = v;
%                     break;
%                 end
%             end
% 
%             if isempty(dataStruct)
%                 error('No valid structure found in file. Check that it contains fields like "FFT" or "LSQ".');
%             end


            % Match the substructure based on the method
            switch methodLower
                case 'fft'
                    substructName = 'FFT';
                case 'full fft'
                    substructName = 'FullFFT';
                case 'lsq'
                    substructName = 'LSQ';
                case 'full lsq'
                    substructName = 'FullLSQ';
                otherwise
                    error('Unexpected method: %s', methodLower);
            end

%             % Validate the structure field
%             if ~isfield(results, substructName)
%                 error('Error: The results structure does not contain "%s" field.', substructName);
%             end

            thisResult = [];
            % If the file contains a top-level variable with same name (like FFT)
            if isfield(temp, substructName)
                thisResult = temp.(substructName);
            else
                % If the file contains wrapper 'results'
                if isfield(temp, 'results') && isstruct(temp.results) && isfield(temp.results, substructName)
                    thisResult = temp.results.(substructName);
                else
                    % Search for any top-level struct that contains DoLP & AoLP
                    fns = fieldnames(temp);
                    for k = 1:numel(fns)
                        candidate = temp.(fns{k});
                        if isstruct(candidate) && isfield(candidate, 'DoLP') && isfield(candidate, 'AoLP')
                            % If the candidate also has the requested substructName as a field (rare), prefer that
                            if isfield(candidate, substructName)
                                thisResult = candidate.(substructName);
                            else
                                thisResult = candidate;
                            end
                            break;
                        end
                    end
                end
            end

            if isempty(thisResult)
                % Auxiliary diagnostic message listing variables found
                varsInFile = strjoin(fieldnames(temp), ', ');
                error(['No usable results structure found in file. Expected a variable named "%s", ' ...
                       'or a wrapper "results.%s", or a struct containing DoLP/AoLP. Variables found: %s'], ...
                       substructName, substructName, varsInFile);
            end


            % Extract DoLP/AoLP
            if isfield(thisResult, 'DoLP') && isfield(thisResult, 'AoLP')
                output1 = thisResult.DoLP;
                output2 = thisResult.AoLP;
            else
                error('The selected file does not contain fields "DoLP" and/or "AoLP".');
            end


        %% in case the choosen method isn't recognized
        otherwise
            error('Error: Unknown HSI method: %s', method);
    end

end
