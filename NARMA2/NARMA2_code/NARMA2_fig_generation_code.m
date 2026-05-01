%% -----------------------------------------------------------------------
%% NARMA-2 Figure: target vs. stochastic estimate (Nres=100, window=1)
%% -----------------------------------------------------------------------

%% 1. Set parameters to match the desired combination
movmean_window  = 1;      % no filtering — best NARMA-2 result
Nres_init       = 100;
Nres            = Nres_init * memorywindowlength;

%% 2. Re-apply filter (trivially no-op at window=1) and re-sample
Active_rec_fin_plot = movmean(ReceptorActives{1,1} / N, [movmean_window 0]);

%% 3. Rebuild test state matrix for Nres=100, window=1
time_index_smol_plot   = START_TIME : TIME_STEP * SAMPLING_PERIOD : STOP_TIME;
steps_per_symbol_smol  = T / TIME_STEP;
sampled_test_plot      = zeros(Nres, num_test_points);
del_idx                = [];

for i = washout1 + num_train_points + washout2 + 1 : num_tot_points
    t_start = (i - 1) * T + offset;
    idx0    = find(time_index_smol_plot >= t_start, 1);
    sidx    = idx0 + (0:Nres-1) * (steps_per_symbol_smol / Nres);
    sidx(sidx > length(time_index_smol_plot)) = [];
    col     = i - (washout1 + num_train_points + washout2);
    if length(sidx) >= Nres
        sampled_test_plot(:, col) = Active_rec_fin_plot(sidx(1:Nres));
    else
        del_idx = [del_idx col];
    end
end
sampled_test_plot(:, del_idx) = [];

%% 4. Re-train readout on same combination (train matrix)
sampled_train_plot = zeros(Nres, num_train_points);
del_idx_tr = [];
for i = washout1 + 1 : washout1 + num_train_points
    t_start = (i - 1) * T + offset;
    idx0    = find(time_index_smol_plot >= t_start, 1);
    sidx    = idx0 + (0:Nres-1) * (steps_per_symbol_smol / Nres);
    sidx(sidx > length(time_index_smol_plot)) = [];
    col     = i - washout1;
    if length(sidx) >= Nres
        sampled_train_plot(:, col) = Active_rec_fin_plot(sidx(1:Nres));
    else
        del_idx_tr = [del_idx_tr col];
    end
end
sampled_train_plot(:, del_idx_tr) = [];

X_tr   = [sampled_train_plot; ones(1, size(sampled_train_plot,2))]';
q_tr   = q(washout1+2 : washout1+num_train_points+1)';
W_plot = pinv(X_tr'*X_tr + 1*eye(size(X_tr,2))) * (X_tr' * q_tr(:));

X_te      = [sampled_test_plot; ones(1, size(sampled_test_plot,2))]';
q_te      = q(washout1+num_train_points+washout2+2 : num_tot_points+1)';
q_te_hat  = X_te * W_plot;

%% 5. Compute NRMSE
nrmse_plot = sqrt(mean((q_te_hat(1:end-2) - q_te(1:end-2)).^2)) / ...
             std(q_te(1:end-2));
fprintf('NARMA-2 NRMSE (Nres=100, W=1): %.4f\n', nrmse_plot);

%% 6. Plot — match style
N_show = 300;   % number of test symbols to display
figure('Position', [100 100 900 320]);
plot(q_te(1:N_show), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
hold on;
plot(q_te_hat(1:N_show), 'r-', 'LineWidth', 1.5);
xlabel('Symbol index');
ylabel('Target and estimated values');
legend('Target', 'MC-PRC (Stochastic)', 'Location', 'best');
title(sprintf('NARMA-2 — Stochastic prediction  |  NRMSE = %.4f', nrmse_plot));
grid on;
xlim([1 N_show]);
set(gca, 'FontSize', 11);