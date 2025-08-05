
function plot_polarization_parameter(wavelengths, data, idxMax, idxMin, max, min, yLabel, titleText)
    fig = figure(); hold on;
    plot(wavelengths, data, 'b--', 'DisplayName','WG Polarizer');
    scatter([wavelengths(idxMax), wavelengths(idxMin)], [max, min], 'b', 'DisplayName', 'Max/Min WG Polarizer');
    text(wavelengths(idxMax), max, sprintf('Max: %.4f', max), 'Color', 'b');
    text(wavelengths(idxMin), min, sprintf('Min: %.4f', min), 'Color', 'b');

    % Labels, title, and formatting
    xlabel('\lambda (nm)');
    ylabel(yLabel);
    title(titleText);
    legend('Location', 'Best');
    grid on; grid minor;
    hold off;

    export_figure(fig, [strrep(titleText, ' ', '_')], 'figures');
end