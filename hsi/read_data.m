function [Data, White, Dark, wavelengths] = read_data(fpath, fname)
    [Data, info] = enviread(strcat(fpath,fname,'.raw'), strcat(fpath,fname,'.hdr')); % Read the hyperspectral data and associated information
    [White, ~] = enviread(strcat(fpath,'WHITEREF_',fname,'.raw'), strcat(fpath,'WHITEREF_',fname,'.hdr')); % Read the white reference data
    [Dark, ~] = enviread(strcat(fpath,'DARKREF_',fname,'.raw'), strcat(fpath,'DARKREF_',fname,'.hdr')); % Read the dark reference data
    wavelengths = str2num(info.Wavelength(2:end-1)); % Extract wavelength information and convert it from string to numeric
end
