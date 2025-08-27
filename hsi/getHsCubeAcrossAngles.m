function cube = getHsCubeAcrossAngles(Wavelength_reflectances) % Multiple hs cubes cause MATLAB memory to run out
        fields = fieldnames(Wavelength_reflectances); % Get all fieldnames (angles)
        % Initialize accumulator with first cube (angle_0)
        firstCube = Wavelength_reflectances.(fields{1});
        cubeSum = zeros(size(firstCube), 'single');
        for f = 1:numel(fields) % Sum across all angles
            cubeSum = cubeSum + Wavelength_reflectances.(fields{f});
        end
        cube = cubeSum ./ numel(fields); % Mean across angles
end