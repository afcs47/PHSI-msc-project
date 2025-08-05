function [sid_out, sam_out, sidsam_out, jmsam_out, ns3_out] = compute_spectral_metrics(reflectance1,reflectance2)

    % Calculate some spectral measures (additional (Add on) hyperspectral toolbox for image processing toolbox necessary)
    % https://www.mathworks.com/help/images/hyperspectral-image-processing.html
    
    % Calculates various spectral similarity metrics
    sid_out = sid(reflectance1,reflectance2); % Spectral Information Divergence
    sam_out = sam(reflectance1,reflectance2); % Spectral Angle Mapper
    sidsam_out = sidsam(reflectance1,reflectance2); % Combined SID and SAM metric
    jmsam_out = jmsam(reflectance1,reflectance2); % Jeffries Matusita-Spectral Angle Mapper 
    ns3_out = ns3(reflectance1,reflectance2); % Normalized Spectral Similarity Score
end
