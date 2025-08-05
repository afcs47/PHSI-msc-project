function set_aolp_colorbar(ax) % 
    caxis(ax, [0 pi]);
    cb = colorbar(ax);
    cb.Ticks = linspace(0, pi, 10);
    cb.TickLabels = arrayfun(@(x) sprintf('%d°', round(rad2deg(x))), cb.Ticks, 'UniformOutput', false);
end
