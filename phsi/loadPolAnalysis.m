function [DoLP_POL, AoLP_POL] = loadPolAnalysis(method, datasetName)  % process POL data according to selected method ( 'Standard', '2nd Order Fourier', '4th Order Fourier')
    
    polDataFile = fullfile(['pol results+figures\' datasetName 'jun.raw'], sprintf('%sjun.raw_polarization_data.mat', datasetName));
    if ~isfile(polDataFile); error('Error: Polarization data file for %s not found! Please run the respective script first or check directories.', datasetName); end 
    
    %% Branch by method 
    switch lower(method)
        case 'standard'
            load(polDataFile, 'DoLP', 'AoLP');
            DoLP_POL = DoLP;
            AoLP_POL = AoLP;

        case '2nd order fourier'
            load(polDataFile, 'DoLP_2nd', 'AoLP_2nd');
            DoLP_POL = DoLP_2nd;
            AoLP_POL = AoLP_2nd;

        case '4th order fourier'
            load(polDataFile, 'DoLP_4th', 'AoLP_4th');
            DoLP_POL = DoLP_4th;
            AoLP_POL = AoLP_4th;
        
        otherwise
            error('Error: Unknown POL method: %s', method);
    end
    
end
