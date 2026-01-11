% Polarization Peanuts Experiment - Expected values

% Clear environment
clear; clc; close all;

% Initialize storage
results = table();

while true
    %% SECTION 1: Setup Parameters
    % Define analyzer angle alpha from 0° to 360° in 10° steps
    alpha_deg = 0:10:360;
    alpha_rad = deg2rad(alpha_deg); % convert to radians

     % Define possible polarizer angles
    pol_angles = [0, 90, 180];
    pol_labels = {'0° (HPL)', '90° (VPL)', '180° (HPL)'};


    % Menu for polarizer angle selection
    [selection1, ok1] = listdlg('PromptString','Select Polarizer Angle:', 'SelectionMode','single', 'ListString',pol_labels);

    if ~ok1
        disp('No selection made. Exiting...');
        break; % exit the loop
    end

    beta = pol_angles(selection1);  % selected polarizer angle in degrees (0 or 180 for HPL; 90 for VPL)
    fprintf('\nSelected polarizer angle: %d°\n', beta);

    while true

        % Define possible QWP angles
        qwp_angles = [0, 30, 45, 60, 90, -30, -45, -60, 180];
        qwp_labels = {'0° (Linear)', '30° (Elliptical)', '45° (Circular)', '60° (Elliptical)', '90° (Linear)', '-30° (Elliptical)', '-45° (Circular)', '-60° (Elliptical)', '180° (Linear)'};
    
        % Menu for QWP angle selection
        [selection2, ok2] = listdlg('PromptString','Select QWP Angle:', 'SelectionMode','single', 'ListString',qwp_labels);
    
        if ~ok2
            disp('No selection made. Exiting...');
            break; % exit the loop
        end
    
        theta = qwp_angles(selection2);  % selected QWP angle in degrees
    
        fprintf('\nSelected QWP angle: %d°\n\n', theta);

    
        %% SECTION 2: Simulate the intensity profile using Jones calculus
    
        % Jones matrix definitions
        jones_polarizer = @(phi) [cosd(phi)^2, cosd(phi)*sind(phi); cosd(phi)*sind(phi), sind(phi)^2];
    
        %jones_qwp = @(phi) [cosd(phi)^2 + 1i*sind(phi)^2, (1 - 1i)*cosd(phi)*sind(phi); (1 - 1i)*cosd(phi)*sind(phi), sind(phi)^2 + 1i*cosd(phi)^2];
        jones_qwp = @(phi) [cos(pi/4) + 1i*sin(pi/4)*cosd(2*phi), 1i*sin(pi/4)*sind(2*phi); 1i*sin(pi/4)*sind(2*phi), cos(pi/4) - 1i*sin(pi/4)*cosd(2*phi)];
    
        % Input light
        E_in = [1; 0];  % [1; 0] for HPL or [0; 1] for VPL

        % Output light
        E_pol = jones_polarizer(beta) * E_in; % light after passing through the polarizer
        E_out = jones_qwp(theta) * E_pol; % light after passing through the QWP 
    
    
        %% SECTION 3: Compute Intensity vs Analyzer Angle
        
        I = zeros(size(alpha_rad));  % initialize intensity array
        for k = 1:length(alpha_rad)
            analyzer = jones_polarizer(alpha_deg(k)); % analyzer at angle
            E_detected = analyzer * E_out;  % detected field
            I(k) = abs(E_detected(1))^2 + abs(E_detected(2))^2; % total intensity
        end
        I = I / max(I); % normalize intensity
    
 %% SECTION 4: Fourier Analysis
        
        % Fit a sinusoid: a0 + a2*cos(2x) + b2*sin(2x)
        fitfun = @(b, x) b(1) + b(2)*cosd(2*x) + b(3)*sind(2*x);
        beta0 = [0.5, 0.5, 0.5]; % initialize fit
        bfit = lsqcurvefit(fitfun, beta0, alpha_deg, I);  % least-squares fit
        
        a0 = bfit(1);
        a2 = bfit(2);
        b2 = bfit(3);
        
        % Compute Stokes parameters from Fourier coefficients
        S0 = 2*a0;
        S1 = 2*a2;
        S2 = 2*b2;
        S3 = NaN;  % cannot be determined from this setup (no info about circular polarization)
    
        % Compute tilt and ellipticity angles
        psi = 0.5 * atan2(S2, S1);
        chi = 0.5 * asin(S3 / S0);
    
        % Display results
        fprintf('Fourier Coefficients:\n');
        fprintf('a0 = %.4f\na1 = %.4f\nb1 = %.4f\n', a0, a2, b2);
        fprintf('\nStokes Parameters:\n');
        fprintf('S0 = %.4f\nS1 = %.4f\nS2 = %.4f\nS3 = %.4f\n', S0, S1, S2, S3);
        fprintf('\nStokes Parameters (normalized):\n');
        fprintf('S0 = %.4f\nS1 = %.4f\nS2 = %.4f\nS3 = %.4f\n', S0/S0, S1/S0, S2/S0, S3/S0);
        fprintf('\nTilt angle: \x03A8 = %.2f deg\n', rad2deg(psi));
        fprintf('Ellipticity angle: \x03c7 = %.2f deg\n', rad2deg(chi));

        %% SECTION 5: Plot Results

        %figure;
        f1 = figure('Name', 'Polarization Analysis', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    
        % Plot intensity vs analyzer angle with fit
        subplot(2,2,1);
        plot(alpha_deg, I, 'o', 'LineWidth', 2);
        xlabel('Analyzer Angle (°)'); ylabel('Normalized Intensity');
        title(sprintf('Intensity vs Analyzer Angle (QWP = %d°)', theta));
        grid on;
        hold on;

        plot(0:1:360, fitfun(bfit, 0:1:360), 'r--', 'LineWidth', 1);
        legend('Data', 'LSQ Fit');
    
        % Polar plot of intensity
        subplot(2,2,2);
        polarplot(alpha_rad, I, 'b-', 'LineWidth', 2);
        title('Polar Plot');
    
        % Polarization ellipse plot
        subplot(2,2,[3 4]);
        X = I .* cos(alpha_rad);
        Y = I .* sin(alpha_rad);
        plot(X, Y, 'm-', 'LineWidth', 2); axis equal; grid on;
        xlabel('I cos(\alpha)'); ylabel('I sin(\alpha)');
        title('Polarization Ellipse');
    
        % Allow user to draw major axis (a)
        fprintf('\nClick two points for the major axis (a)...\n');
        [xa, ya] = ginput(2); % select points
        hold on; plot(xa, ya, 'r-', 'LineWidth', 2);
        dx_a = xa(2) - xa(1); dy_a = ya(2) - ya(1); % distance between choosen points
        a = sqrt(dx_a^2 + dy_a^2);  % length of major axis
        psi_manual = atan2(dy_a, dx_a);  % angle of major axis with horizontal axis; ]0; pi] -> atan2 [–π, π]
    
        % Allow user to draw minor axis (b)
        fprintf('Click two points for the minor axis (b)...\n');
        [xb, yb] = ginput(2); % select points
        plot(xb, yb, 'g-', 'LineWidth', 2);
        dx_b = xb(2) - xb(1); dy_b = yb(2) - yb(1); % distance between choosen points
        b = sqrt(dx_b^2 + dy_b^2); % length of minor axis
    
        chi_manual = atan(sqrt(b / a)); % ellipticity angle from axes ratio; ]-pi/4; pi/4] -> atan [–π/2, π/2]
    
        % Display manual measurement results
        fprintf('\nManual Measurement Results:\n');
        fprintf('a = %.4f, b = %.4f\n', a, b);
        fprintf('Tilt angle \x03A8 = %.2f°\n', rad2deg(psi_manual));
        fprintf('Ellipticity angle \x03c7 = %.2f°\n\n', rad2deg(chi_manual));
    
        
        %S0 = 1; % normalized
        S3_manual = S0 * sin(2*chi_manual);
        
        % Warn if nearly circular
        if abs(a - b)/max(a,b) < 0.05 % 5% tolerance
            fprintf('Polarization nearly circular, S3 sign undetermined\n');
            S3_manual = abs(S3_manual); % only magnitude
        end
        fprintf('S3_manual = %.2f\n', S3_manual);
        fprintf('S3_manual (normalized) = %.2f\n\n', S3_manual/S0);
    
        % Collect results
        results = [results; table(beta, theta, a0, a2, b2, S0, S1, S2, S3,...
            rad2deg(psi), rad2deg(chi), ...
            rad2deg(psi_manual), rad2deg(chi_manual), S3_manual,...
            'VariableNames', {'Beta (°)', 'Theta (°)', 'a0', 'a2', 'b2', 'S0', 'S1', 'S2', 'S3', ...
            'Psi (°)', 'Chi (°)',...
            'Psi_manual  (°)', 'Chi_manual (°)', 'S3_manual'})];

        % Pause and return to menu
        pause(0.5);
 
    end

end
% Save results
writetable(results, fullfile('E:\5o ano\Tese\pol_setup_peanuts', ['Polarization_Analysis_Expected_Values_', char(datetime('now','Format','ddMMyyyy''_''HHmmss')),'.xlsx']));
disp('Results saved.');