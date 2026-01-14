% Polarization Peanuts Experiment - Experimental values + fourier analysis

% Clear environment
clear; clc; close all;

%main
[raw_data, num_sets, I_background, fullpath] = loadData();
[pathname,~,~] = fileparts(fullpath); % get directory


% Load theoretical results
[theory_file, theory_path] = uigetfile({'*.xlsx;*.xls'}, 'Select the Excel file with expected results');
if isequal(theory_file,0)
    warning('No file selected. Continuing without expected results for comparison.');
    theory_data = table();
else
    theory_fullpath = fullfile(theory_path, theory_file);
    theory_data = readtable(theory_fullpath);
end

% Initialize storage
results = table();

for i = 1:num_sets
    try
        set_idx =i;
        results = analyzeSet(i, raw_data, I_background, results, theory_data);
       
        pause(0.5);
    catch ME
        warning('Error analyzing set %d: %s', i, ME.message);
        continue; % Save results if error occurs
    end    
end

% Save results
%writetable(results, fullfile(pathname, ['Polarization_Analysis_Results_', char(datetime('now','Format','ddMMyyyy''_''HHmmss')),'.xlsx']));

output_file = fullfile(pathname, ['Polarization_Analysis_Results_', char(datetime('now','Format','ddMMyyyy''_''HHmmss')), '.xlsx']);
writetable(results, output_file, 'Sheet', 'Experimental');

if ~isempty(theory_data)
    writetable(theory_data, output_file, 'Sheet', 'Theoretical');
end

disp('Results saved.');




%% Functions
function [raw_data, num_sets, I_background, fullpath] = loadData()
    % Prompt user to select Excel file
    [filename, pathname] = uigetfile({'*.xlsx;*.xls'}, 'Select the Excel file with the experimental data');
    if isequal(filename, 0)
        error('No file selected. Exiting...');
    end
    fullpath = fullfile(pathname, filename);
    
    raw_data = readmatrix(fullpath,'Sheet',3, 'NumHeaderLines', 4); % skip headers %sheet 2 for 05/05 file and sheet 3 for 29/05 file

    % Check number of measurement sets in columns (8 columns per set for 29/05 file and 6 for 05/05 file)
    num_sets = size(~isnan(raw_data), 2) / 8; % excludes NaN values from the rest of the sheet
    %num_sets = 80 / 8;
    
    I_background = mean([0.128, 0.358]); % values taken after covering the detector ~ 0.243 mV

end

function [theta, alpha, I, dI] = extractMeasurementSet(raw_data, set_idx, I_background)
    x=8; %x = 6 for 05/05 file and 8 for 29/05
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


function [bfit, stderr] = fitFourier(alpha, I, dI)
    fitfun = @(b,x) b(1) + b(2)*cosd(2*x) + b(3)*sind(2*x); % fitted sinusoidal function
    beta0 = [0.5, 0.5, 0.5]; % initial guess
    if isempty(dI)
        bfit = lsqcurvefit(fitfun, beta0, alpha, I);
        stderr = [];
        %plot(0:1:360, fitfun(bfit, 0:1:360), 'g--', 'LineWidth', 1.5, 'DisplayName', 'lsqnonlin fit');
    else
        weighted_residuals = @(params) (fitfun(params, alpha) - I) ./ dI; % Weighted residual function
        options = optimoptions('lsqnonlin', 'Display', 'off'); % Create options for the nonlinear least squares solver; 'Display', 'off' disables output during optimization; 
                                                               % 'Algorithm','levenberg-marquardt' or 'Algorithm','trust-region-reflective' have the same results
        [bfit, ~, ~, ~, ~, ~, J] = lsqnonlin(weighted_residuals, beta0, [], [], options); % Perform nonlinear least squares fitting using weighted residuals;
                                                                                                % No lower or upper bounds are specified ([]), options are passed in;
                                                                                                % J is the Jacobian matrix of partial derivatives at the solution
        ci = nlparci(bfit, weighted_residuals(bfit), 'jacobian', J); % Returns an estimate of the 95% confidence intervals (CI) for the fitted nonlinear least-squares parameter beta using the Jacobian and the residuals evaluated at the solution
        stderr = (ci(:,2) - ci(:,1)) / (2*1.96); % approximate 1σ from 95% CI; 1.96 is the z-score for a 95% confidence interval in a normal distribution
        
        %plot(0:1:360, fitfun(bfit, 0:1:360), 'r--', 'LineWidth', 1, 'DisplayName', 'lsqnonlin fit w/ error');
    end
end

% function [bfit, gof] = fitFourierDiagnostic(alpha,I,dI) %check if a1 and b1 are ~0 -> fit too far off of the actual points for the verification to make sense
%     fitfunc = fittype('a0 + a1*cosd(x) + b1*sind(x) + a2*cosd(2*x) + b2*sind(2*x)');
%     opts = fitoptions('Method','NonlinearLeastSquares');
%     opts.StartPoint = [0.5 0 0 0.5 0];
%     opts.Weights = 1./(dI.^2);
% 
%     [f,gof] = fit(alpha(:),I(:),fitfunc,opts);
%     bfit = coeffvalues(f);
% end

function [S, dS] = computeStokes(bfit, stderr)
    S0 = 2 * bfit(1); S1 = 2 * bfit(2); S2 = 2 * bfit(3); S3 = NaN;
    dS0 = 2 * stderr(1); dS1 = 2 * stderr(2); dS2 = 2 * stderr(3); dS3 = NaN;
    S = [S0, S1, S2, S3];
    dS = [dS0, dS1, dS2, dS3];
end

function [psi, dpsi, chi, dchi] = computeAngles(S, dS)
%     S0 = S(1); S1 = ; S2 = S(3); S3 = S(4);
%     dS0 = dS(1); dS1 = dS(2); dS2 = dS(3); dS3 = dS(4);

    psi = 0.5 * atan2(S(2+1), S(1+1));
    dpsi = 0.5 * sqrt((S(1+1)^2 * dS(2+1)^2 + S(2+1)^2 * dS(1+1)^2) / (S(1+1)^2 + S(2+1)^2)^2);

    chi = 0.5 * asin(S(3+1) / S(0+1)); % = NaN in this setup
    dchi = NaN; % because S3 = NaN
end

function [DoLP, dDoLP, AoLP, dAoLP] = computePolMetrics(S, dS)
%     S0 = S(1); S1 = S(2); S2 = S(3);
%     dS0 = dS(1); dS1 = dS(2); dS2 = dS(3);

    DoLP = sqrt(S(1+1)^2 + S(2+1)^2) / S(0+1);
    dDoLP = sqrt((S(1+1)^2 * dS(1+1)^2 + S(2+1)^2 * dS(2+1)^2)/(S(1+1)^2 + S(2+1)^2) + (S(1+1)^2 + S(2+1)^2) * dS(0+1)^2 / S(0+1)^2) / S(0+1);

    AoLP = rad2deg(0.5 * atan2(S(2+1), S(1+1)));
    dAoLP = rad2deg(0.5 * sqrt((S(1+1)^2 * dS(2+1)^2 + S(2+1)^2 * dS(1+1)^2) / (S(1+1)^2 + S(2+1)^2)^2));
end

function [a0, a2, b2, a4, b4] = extractFourierFromFFT(I)
    N = length(I);
    I_fft = fft(I) / N;
    a0 = real(I_fft(1));
    a2 = 2 * real(I_fft(3));
    b2 = -2 * imag(I_fft(3));
    a4 = 2*real(I_fft(5)); b4 = -2*imag(I_fft(5)); %for diagnostics
end


function [psi, dpsi, chi, dchi] = manualEllipsePlot(alpha, I, set_idx, theta)
    figure('Name', sprintf('Ellipse Measurement - Set %d (θ=%d°)', set_idx, theta));
    X = I .* cosd(alpha); Y = I .* sind(alpha);
    plot(X, Y, 'k-', 'LineWidth', 1.5); axis equal; grid on;
    title('Select major (red) and minor (green) axes');
    xlabel('I cos(\alpha)'); ylabel('I sin(\alpha)'); hold on;
    
    dx=0.02; %dx=dy, visually estimated value
    
    disp('Click two points for major axis');
    [xa, ya] = ginput(2); plot(xa, ya, 'r-', 'LineWidth', 2);
%     dx_a = xa(2) - xa(1); dy_a = ya(2) - ya(1); % distance between choosen points
%     a = sqrt(dx_a^2 + dy_a^2);  % length of major axis
    a = norm([xa(2)-xa(1), ya(2)-ya(1)]); % length of major axis from the distance between choosen points
    da = sqrt(2)*dx;
    fprintf('Major axis estimate: a=%.4f±%.4f\n',a,da);
    psi = atan2(ya(2) - ya(1), xa(2) - xa(1)); % angle of major axis with horizontal axis; ]0; pi] -> atan2 [–π, π]
    dpsi = sqrt(2)*dx/a;

    disp('Click two points for minor axis');
    [xb, yb] = ginput(2); plot(xb, yb, 'g-', 'LineWidth', 2);
%     dx_b = xb(2) - xb(1); dy_b = yb(2) - yb(1); % distance between choosen points
%     b = sqrt(dx_b^2 + dy_b^2);  % length of minor axis
    b = norm([xb(2)-xb(1), yb(2)-yb(1)]); % length of minor axis from the distance between choosen points
    db = da; %db = sqrt(2)*dy , dy=dx;
    fprintf('Minor axis estimate: b=%.4f±%.4f\n\n',b,db);
    chi = atan(sqrt(b/a)); % ellipticity angle from axes ratio; ]-pi/4; pi/4] -> atan [–π/2, π/2]
    dchi = sqrt(0.5*(b/a + a/b)/(b+a))*da;
end

function [fourier_fit_params, gof, fitted_curve] = fitFourierCurvefit(alpha, I, dI)
    % Ensure inputs are column vectors
    alpha = alpha(:);
    I = I(:);
    dI = dI(:);

    ft = fittype('a0 + a2*cosd(2*x) + b2*sind(2*x) + a4*cosd(4*x) + b4*sind(4*x)');

    % Fit options: weight by 1 / (dI^2)
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.StartPoint = [0.5,0.5, 0.5, 0.1, 0.1]; % initial guess
    opts.Weights = 1 ./ (dI.^2);

    % Perform the fit
    [fitted_curve, gof] = fit(alpha, I, ft, opts); % no ';' to provide Goodness-of-fit statistics

    % Extract coefficients: [a0, a2, b2]
    coeffs = coeffvalues(fitted_curve);
    fourier_fit_params = coeffs;
end

function results = analyzeSet(set_idx, raw_data, I_background, results, theory_data)
    [theta, alpha, I, dI] = extractMeasurementSet(raw_data, set_idx, I_background);

    if (exist('theory_data','var') && ~isempty(theory_data))
        % Match by theta
        if any(strcmp('Beta___', theory_data.Properties.VariableNames))
            idx = theory_data.("Theta___") == theta;
        
        end
        if any(idx)
            theory_row = theory_data(find(idx,1),:);
            %Theory
            a0_expected = theory_row.a0;
            a2_expected = theory_row.a2;
            b2_expected = theory_row.b2;
        end
    end 

    % Fits
    %LSQ
    [bfit, ~] = fitFourier(alpha, I, []);
    [bfit_weighted, stderr] = fitFourier(alpha, I, dI); %2nd order harmonics
    [bfit_curvefit, gof, fitted_curve_obj] = fitFourierCurvefit(alpha, I, dI); %4th order harmonics
    %[fit_full,gof_full] = fitFourierDiagnostic(alpha,I,dI); %verify a1 and b1 are null (ideal conditions) - fit too far off the actual points for it to be considered
    
    %FFT
    [a0,a2,b2,a4,b4] = extractFourierFromFFT(I);
    
    % Polarization parameters
    [S, dS] = computeStokes(bfit_weighted, stderr);
    [psi, dpsi, chi, dchi] = computeAngles(S, dS);
    [DoLP, dDoLP, AoLP, dAoLP] = computePolMetrics(S, dS);

    % Manual ellipse
    [psi_m,dpsi_m,chi_m,dchi_m] = manualEllipsePlot(alpha,I,set_idx, theta);
    S3_m = S(1) * sin(2*chi_m);
    dS3_m = sqrt(( sin(2*chi_m)*dS(1) )^2 + ( 2*S(1)*cos(2*chi_m)*dchi_m )^2);

    % Console summary
    fprintf('\n--- Set %d | θ = %.0f° ---\n',set_idx,theta);
    if exist('a0_expected','var')
        fprintf('Expected values: a0=%.4f, a2=%.4f, b2=%.4f\n', a0_expected, a2_expected, b2_expected);
    end
    fprintf('LSQ 2nd order Fit wo/ errors: a0=%.4f, a2=%.4f, b2=%.4f\n', bfit(1), bfit(2), bfit(3));
    fprintf('LSQ 2nd order Fit: a0=%.4f±%.4f, a2=%.4f±%.4f, b2=%.4f±%.4f\n', bfit_weighted(1), stderr(1), bfit_weighted(2), stderr(2), bfit_weighted(3), stderr(3));
    fprintf('LSQ 4th order Fit: a0=%.4f, a2=%.4f, b2=%.4f, R²=%.4f\n', bfit_curvefit(1), bfit_curvefit(2), bfit_curvefit(3), gof.rsquare);
    %fprintf('LSQ Full Fourier Fit: a0=%.4f, a1=%.4f, b1=%.4f, a2=%.4f, b2=%.4f, R²=%.4f\n', fit_full(1), fit_full(2), fit_full(3),fit_full(4), fit_full(5), gof_full.rsquare);
    fprintf('FFT Coefficients: a0=%.4f, a2=%.4f, b2=%.4f, a4=%.4f, b4=%.4f\n\n', a0, a2, b2, a4, b4);

    fprintf('For LSQ 2nd order Fit coefficients\n');
    
    fprintf('Stokes: S0=%.3f±%.3f, S1=%.3f±%.3f, S2=%.3f±%.3f\n', S(1),dS(1),S(2),dS(2),S(3),dS(3));
    fprintf('\x03A8 (Fourier) = %.2f° ± %.2f°\n',rad2deg(psi),rad2deg(dpsi));
    fprintf('\x03c7 (Fourier) = %.2f° ± %.2f°\n',rad2deg(chi),rad2deg(dchi));
    fprintf('\x03A8 (Manual)  = %.2f° ± %.2f°\n',rad2deg(psi_m),rad2deg(dpsi_m));
    fprintf('\x03c7 (Manual)  = %.2f° ± %.2f°\n',rad2deg(chi_m),rad2deg(dchi_m));
    fprintf('S3 (Manual)  = %.2f° ± %.2f°\n\n',S3_m,dS3_m);

    fprintf('DoLP = %.4f ± %.4f\n', DoLP, dDoLP);
    fprintf('AoLP = %.4f° ± %.4f°\n\n', AoLP, dAoLP);


    % Plot summary
%     figure('Name',sprintf('Comparison of fits for QWP - Set %d (θ=%d°)', set_idx, theta));
%     errorbar(alpha,I,dI,'ko','MarkerFaceColor', 'b', 'DisplayName', 'Data'); hold on;
%     af = 0:1:360; %"alpha"
%     plot(af,bfit(1)+bfit(2)*cosd(2*af)+bfit(3)*sind(2*af),'g--', 'LineWidth', 1.5, 'DisplayName', 'LSQ 2nd order fit wo/ error');
%     plot(af, bfit_weighted(1)+bfit_weighted(2)*cosd(2*af)+bfit_weighted(3)*sind(2*af), 'r--', 'LineWidth', 1, 'DisplayName', 'LSQ 2nd order fit');
%     plot(af, fitted_curve_obj(af), 'b--', 'LineWidth', 1.5, 'DisplayName', 'LSQ 4th order fit');
%     %plot(af,fit_full(1)+fit_full(2)*cosd(af)+fit_full(3)*sind(af)+ fit_full(4)*cosd(2*af)+fit_full(5)*sind(2*af),'g--', 'LineWidth', 1.5, 'DisplayName',  'LSQ Full Fourier fit');
%     plot(af,a0+a2*cosd(2*af)+b2*sind(2*af),'m--', 'LineWidth', 1.5, 'DisplayName', 'FFT 2nd order');
%     plot(af,a0+a2*cosd(2*af)+b2*sind(2*af)+ a4*cosd(4*af)+b4*sind(4*af),'c--', 'LineWidth', 1.5, 'DisplayName', 'FFT 4th order');
%     if (exist('theory_data','var') && ~isempty(theory_data))
%     plot(af, a0_expected + a2_expected*cosd(2*af) + b2_expected*sind(2*af), 'k--', 'LineWidth', 1, 'DisplayName', 'Expected behaviour');
%     end
%     legend;
%     title(sprintf('Comparison of fits for QWP - Set %d (θ=%d°)', set_idx, theta));
%     xlabel('\alpha (°)'); ylabel('Normalized I'); grid on;
%     ylim([-0.05, 1.15]);
%     xlim([0,365]);

    % Plot summary – split LSQ / FFT vs Expected
    figure('Name', sprintf('Comparison of fits for QWP - Set %d (θ=%d°)', set_idx, theta));
    af = 0:1:360; % alpha grid
    %Expected vs LSQ fits
    subplot(2,1,1)
    errorbar(alpha, I, dI, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Data'); hold on;   
    plot(af, bfit(1) + bfit(2)*cosd(2*af) + bfit(3)*sind(2*af), 'g--', 'LineWidth', 1.5, 'DisplayName', 'LSQ 2nd order (unweighted)');
    plot(af, bfit_weighted(1) + bfit_weighted(2)*cosd(2*af) + bfit_weighted(3)*sind(2*af), 'r--', 'LineWidth', 1.2, 'DisplayName', 'LSQ 2nd order (weighted)');  
    plot(af, fitted_curve_obj(af),'b--', 'LineWidth', 1.5, 'DisplayName', 'LSQ 4th order');
    if exist('a0_expected','var')
        plot(af, a0_expected + a2_expected*cosd(2*af) + b2_expected*sind(2*af), 'k--', 'LineWidth', 1, 'DisplayName', 'Expected behaviour');
    end
    title('Expected vs LSQ fits');
    ylabel('Normalized I');
    grid on;
    ylim([-0.05, 1.15]);
    xlim([0, 365]);
    legend('Location','best');
    
    % Expected vs FFT
    subplot(2,1,2)
    errorbar(alpha, I, dI, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Data'); hold on;
    plot(af, a0 + a2*cosd(2*af) + b2*sind(2*af), 'm--', 'LineWidth', 1.5, 'DisplayName', 'FFT 2nd order');
    plot(af, a0 + a2*cosd(2*af) + b2*sind(2*af) + a4*cosd(4*af) + b4*sind(4*af), 'c--', 'LineWidth', 1.5, 'DisplayName', 'FFT 4th order');
    if exist('a0_expected','var')
        plot(af, a0_expected + a2_expected*cosd(2*af) + b2_expected*sind(2*af),'k--', 'LineWidth', 1, 'DisplayName', 'Expected behaviour');
    end
    title('Expected vs FFT reconstructions');
    xlabel('\alpha (°)');
    ylabel('Normalized I');
    grid on;
    ylim([-0.05, 1.15]);
    xlim([0, 365]);
    legend('Location','best');


    %% Collect results
    results = [results; table(theta, a0, a2, b2, a4, b4, S(1),dS(1),S(2),dS(2),S(3),dS(3), ...
        rad2deg(psi),rad2deg(dpsi), rad2deg(chi), rad2deg(dchi), ...
        rad2deg(psi_m), rad2deg(dpsi_m), rad2deg(chi_m), rad2deg(dchi_m), S3_m, dS3_m, ...
        'VariableNames', {'Theta (°)', 'a0', 'a2', 'b2', 'a4', 'b4', 'S0','dS0','S1','dS1','S2','dS2' ...
        'Psi (°)', 'dPsi (°)', 'Chi (°)', 'dChi (°)'...
        'Psi_manual (°)', 'dPsi_m (°)', 'Chi_manual (°)', 'dChi_m (°)', 'S3_manual', 'dS3_m'})];

end

