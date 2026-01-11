% Polarization Peanuts Experiment - Experimental values

% Clear environment
clear; clc; close all;

%main
[raw_data, num_sets, I_background, fullpath] = loadData();
[pathname,~,~] = fileparts(fullpath); % get directory

% Initialize storage
results = table();

for i = 1:num_sets
    try
        results = analyzeSet(i, raw_data, I_background, results);
        pause(0.5);
    catch ME
        warning('Error analyzing set %d: %s', i, ME.message);
        break; % Save results if error occurs
    end    
end

% Save results
writetable(results, fullfile(pathname, ['Polarization_Analysis_Results_', char(datetime('now','Format','ddMMyyyy''_''HHmmss')),'.xlsx']));
disp('Results saved.');




%% Functions
function [raw_data, num_sets, I_background, fullpath] = loadData()
    % Prompt user to select Excel file
    [filename, pathname] = uigetfile({'*.xlsx;*.xls'}, 'Select the Excel file with the experimental data');
    if isequal(filename, 0)
        error('No file selected. Exiting...');
    end
    fullpath = fullfile(pathname, filename);
    
    raw_data = readmatrix(fullpath,'Sheet',2, 'NumHeaderLines', 4); % skip headers %sheet 2 for 05/05 file and sheet 3 for 29/05 file

    % Check number of measurement sets in columns (8 columns per set for 29/05 file and 6 for 05/05 file)
    num_sets = size(~isnan(raw_data), 2) / 6; % excludes NaN values from the rest of the sheet
    %num_sets = 80 / 8;
    
    I_background = mean([0.128, 0.358]); % values taken after covering the detector ~ 0.243 mV

end

function [theta, alpha, I, dI] = extractMeasurementSet(raw_data, set_idx, I_background)
    x=6; %x = 6 for 05/05 file and 8 for 29/05
    col_start = (set_idx - 1) * x + 1;
    % Extract values for the current set
    theta = raw_data(1, col_start);
    d_theta = raw_data(1, col_start + 1);
    alpha = raw_data(:, col_start + 2);
    d_alpha = raw_data(1, col_start + 3);
    I = raw_data(:, col_start + 4);
    dI = raw_data(:, col_start + 5);

    % Remove NaNs or incomplete rows (there are only values on the excel sheet until row 37)
    valid_rows = ~isnan(alpha) & ~isnan(I);
    alpha = alpha(valid_rows);
    I = I(valid_rows) - I_background; % Subtract the background from all intensity values
    I(I < 0) = 0; % Optionally clip negatives to zero
    dI = 0.5*dI(valid_rows)/max(I); % Normalize uncertainties (0.5 added since value in the excel sheet corresponds to 2*u, the total oscilattion of the value presented, instead of the usual uncertainty of the measuring device)
    I = I / max(I); % Normalize intensity

end

function results = analyzeSet(set_idx, raw_data, I_background, results)
    [theta, alpha_deg, I, dI] = extractMeasurementSet(raw_data, set_idx, I_background);
    alpha_rad = deg2rad(alpha_deg);
    beta = 0;
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
    

end


