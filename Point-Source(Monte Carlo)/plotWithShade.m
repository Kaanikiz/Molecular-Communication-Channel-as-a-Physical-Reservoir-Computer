function plotWithShade(t, y, s, labelText)

    t = t(:);
    y = y(:);
    s = s(:);

    valid = isfinite(t) & isfinite(y) & isfinite(s);
    t = t(valid);
    y = y(valid);
    s = s(valid);

    fill([t; flipud(t)], [y-s; flipud(y+s)], ...
        [0.8 0.8 0.8], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.25, ...
        'HandleVisibility', 'off');

    plot(t, y, 'LineWidth', 1.5, 'DisplayName', labelText)
end