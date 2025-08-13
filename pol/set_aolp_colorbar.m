function set_aolp_colorbar(ax, min_angle, max_angle) % 
    caxis(ax, [deg2rad(min_angle) deg2rad(max_angle)]);
    cb = colorbar(ax);
    cb.Ticks = linspace(deg2rad(min_angle), deg2rad(max_angle), 10);
    cb.TickLabels = arrayfun(@(x) sprintf('%d°', round(rad2deg(x))), cb.Ticks, 'UniformOutput', false);
end
