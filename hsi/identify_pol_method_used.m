function method = identify_pol_method_used(name)
method = 'Unknown';
    if contains(name, 'Reflectances')
        method = 'Reflectances';
    elseif contains(name, 'Standard_Stokes')
        method = 'Standard';
    elseif contains(name, 'Fourier_LSQ')
        method = 'Fourier_LSQ';
    elseif contains(name, 'Fourier_Full_LSQ')
        method = 'Fourier_Full_LSQ';
    elseif contains(name, 'Fourier_FFT')
        method = 'Fourier_FFT';
    elseif contains(name, 'Fourier_Full_FFT')
        method = 'Fourier_Full_FFT';
    elseif contains(name, 'SPIE_based')
        method = 'SPIE';
    elseif contains(name, 'SPIE_Simplified')
        method = 'SPIE_Simple';
    end
end