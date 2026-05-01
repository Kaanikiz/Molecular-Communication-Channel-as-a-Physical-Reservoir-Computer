load('pointApproxValidation_results.mat')

%% Plot-only snippet with corr/NRMSE textbox only
% Assumes these variables already exist:
% t_det, b_det, t_point, B_point, finiteResults

%% Point-TX / finite-RX mean and metrics

b_point_mean = mean(B_point, 2, 'omitnan');
b_point_std  = std(B_point, 0, 2, 'omitnan');

% Interpolate analytical trace onto Smoldyn time grid
b_det_point = interp1(t_det, b_det, t_point, 'linear', 'extrap');

metrics_point = computeCorrNRMSE_local(b_det_point, b_point_mean);

%% Finite-TX / finite-RX mean and metrics

for rr = 1:length(finiteResults)

    if ~isfield(finiteResults(rr), 'mean') || isempty(finiteResults(rr).mean)
        finiteResults(rr).mean = mean(finiteResults(rr).B, 2, 'omitnan');
    end

    if ~isfield(finiteResults(rr), 'std') || isempty(finiteResults(rr).std)
        finiteResults(rr).std = std(finiteResults(rr).B, 0, 2, 'omitnan');
    end

    b_det_fin = interp1(t_det, b_det, finiteResults(rr).t, 'linear', 'extrap');

    finiteResults(rr).metrics = computeCorrNRMSE_local( ...
        b_det_fin, finiteResults(rr).mean);
end

%% Plot

figure('Color','w');
hold on; grid on; box on;

% Colors
col_analytical = [0 0 0];             % black
col_pointTX    = [0 0.4470 0.7410];   % blue
col_finiteTX   = [0.8500 0.3250 0.0980]; % red/orange

% Analytical deterministic trace
plot(t_det, b_det, '-', ...
    'Color', col_analytical, ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Analytical point model');

% Point-TX / finite-RX Smoldyn trace: blue line + blue shaded std
plotWithShade_local(t_point, b_point_mean, b_point_std, ...
    col_pointTX, ...
    'Smoldyn point-TX / finite-RX mean');

% Finite-TX / finite-RX Smoldyn trace: red line + red shaded std
for rr = 1:length(finiteResults)
    plotWithShade_local(finiteResults(rr).t, finiteResults(rr).mean, finiteResults(rr).std, ...
        col_finiteTX, ...
        sprintf('Smoldyn finite-TX r_{tx}=%.1f \\mum / finite-RX mean', finiteResults(rr).r_tx));
end

xlabel('Time (s)', 'FontSize', 13);
ylabel('Receptor occupancy ratio b(t)', 'FontSize', 13);
title('Point approximation validation with 20-pulse input', ...
    'FontSize', 14, 'FontWeight', 'bold');

legend('Location','southeast', 'FontSize', 10);

%% Textbox

txt = sprintf(['Analytical vs. Smoldyn point-TX / finite-RX:\n' ...
               '  corr = %.4f\n' ...
               '  NRMSE = %.4f\n'], ...
               metrics_point.corr, metrics_point.nrmse);

for rr = 1:length(finiteResults)
    txt = sprintf(['%s\nAnalytical vs. Smoldyn finite-TX / finite-RX:\n' ...
                   '  corr = %.4f\n' ...
                   '  NRMSE = %.4f\n'], ...
                   txt, ...
                   finiteResults(rr).metrics.corr, ...
                   finiteResults(rr).metrics.nrmse);
end

annotation('textbox', [0.50 0.23 0.50 0.40], ...
    'String', txt, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 10, ...
    'Interpreter', 'tex');

xlim([0 max(t_det)]);
ylim padded;

savefig('pointApproxValidation_corr_nrmse_colored.fig');
saveas(gcf, 'pointApproxValidation_corr_nrmse_colored.png');

%% Helper functions

function metrics = computeCorrNRMSE_local(b_ref, b_test)

    b_ref = b_ref(:);
    b_test = b_test(:);

    valid = isfinite(b_ref) & isfinite(b_test);
    b_ref = b_ref(valid);
    b_test = b_test(valid);

    b_ref_centered = b_ref - mean(b_ref);
    b_test_centered = b_test - mean(b_test);

    metrics.corr = sum(b_ref_centered .* b_test_centered) / ...
        sqrt(sum(b_ref_centered.^2) * sum(b_test_centered.^2));

    metrics.nrmse = sqrt(mean((b_ref - b_test).^2)) / ...
        (max(b_ref) - min(b_ref));
end

function plotWithShade_local(t, y, s, lineColor, labelText)

    t = t(:);
    y = y(:);
    s = s(:);

    valid = isfinite(t) & isfinite(y) & isfinite(s);
    t = t(valid);
    y = y(valid);
    s = s(valid);

    % Same-color uncertainty band
    fill([t; flipud(t)], [y-s; flipud(y+s)], ...
        lineColor, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.28, ...  
        'HandleVisibility', 'off');

    % Mean trace
    plot(t, y, '-', ...
        'Color', lineColor, ...
        'LineWidth', 2.0, ...
        'DisplayName', labelText);
end