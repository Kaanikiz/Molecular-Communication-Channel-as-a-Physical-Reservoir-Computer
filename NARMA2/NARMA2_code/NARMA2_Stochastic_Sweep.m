%% NARMA2_Stochastic_Sweep.m
%
% Stochastic parameter sweep for the NARMA2 benchmark task, mirroring the
% structure of NARMA_Stochastic_Sweep.m (NARMA10 version) exactly.
%
% NARMA2 recurrence (order-2 NARMA):
%   q(n+1) = 0.4*q(n) + 0.4*q(n)*q(n-1) + 0.6*u(n)^3 + 0.1
%
% Key differences from the NARMA10 sweep:
%   - NARMA order: 2 instead of 10.
%   - Warmup in narma_len: +2 instead of +10.
%   - Generation loop starts at n = 2 instead of n = 10.
%%

%% Initialize the MATLAB environment
close all
clear

rng(0); % Fix random seed for reproducibility (numerical model and Smoldyn)

smold           = 1;               % 1 = include Smoldyn path; 0 = numerical only
smold_only_data = 1;               % 1 = read existing Smoldyn output (no new run)
datadirname     = 'NARMA2_5000_T1_koff1'; % Directory with pre-existing Smoldyn data
                                          % (only used when smold=1 & smold_only_data=1)

movmean_on                 = 1;    % Apply causal moving-average to Smoldyn signal
movmean_window             = 500;  % Moving-average window (samples)
memorywindowlength         = 1;    % Number of symbol intervals to look back (virtual nodes)
movmean_readout_test       = 0;    % Apply moving average to readout predictions?
movmean_readout_test_window = 5;   % Window for readout moving average

plot_in_loop = false;  % true = show intermediate figures inside the sweep loop


%% Primary & Secondary Swept Parameters
% -------------------------------------------------------------------------
% param1 is swept in the inner loop; param2 in the outer loop.
% -------------------------------------------------------------------------
param1_name   = 'movmean_window';
param1_values = 1;

param2_name   = 'Nres_init';
param2_values = 100;


%% 1. Base Parameter Setup
% Receptor / binding parameters
N     = 500;       % Total receptors at receiver
k_on  = 1e-18;     % Binding rate constant [m^3/(mol·s)]
k_off = 1;         % Unbinding rate constant [1/s]
KD    = k_off / k_on; % Dissociation constant

% Time parameters
T = 1;             % Symbol duration [s]

% Data-split lengths
washout1         = 250;   % Washout symbols before training
num_train_points = 2000;  % Training symbols
wheretostarttest = 2500;
washout2         = wheretostarttest - (washout1 + num_train_points);
num_test_points  = 2000;  % Test symbols

num_tot_points = washout1 + washout2 + num_train_points + num_test_points;

% Channel parameters
distance = 10e-6;   % Tx–Rx distance [m]
D        = 1e-11;   % Diffusion coefficient [m^2/s]

% Input normalization
N_min = 100;
N_max = 3000;

% Regularization
lambda      = 1e-10; % Ridge regularization (deterministic model)
lambda_smol = 1;     % Ridge regularization (Smoldyn model)

% Reservoir size
Nres_init = 100;
Nres      = Nres_init * memorywindowlength;

% Numerical integration step
dt = 0.001;

% Channel memory window
memlengthsweep = 100;
Tpeak          = distance^2 / (6 * D);
memory_length  = round((memlengthsweep * Tpeak) / T) * T;
offset         = 0;


%% Pre-allocate 2-D result arrays
nP1 = length(param1_values);
nP2 = length(param2_values);

nrmse_train_NUM_array  = zeros(nP1, nP2);
nrmse_test_NUM_array   = zeros(nP1, nP2);
nrmse_train_SMOL_array = zeros(nP1, nP2);
nrmse_test_SMOL_array  = zeros(nP1, nP2);

% Automatic log-scale detection for heatmap axes
minP1   = min(param1_values); maxP1 = max(param1_values);
useLogX = (minP1 > 0) && (maxP1 / minP1 > 50);
minP2   = min(param2_values); maxP2 = max(param2_values);
useLogY = (minP2 > 0) && (maxP2 / minP2 > 50);


%% 2. Generate the NARMA2 Time Series (Once, Before Sweep)
tic
% NARMA2 lag = 2 → warmup = 2 samples
narma_order = 2;
narma_len   = num_tot_points + narma_order + 1; % +2 warmup, +1 for next-step targets

% Random i.i.d. input in [0, 0.5]
u = 0.5 * rand(narma_len, 1);

% Pre-allocate output (zeros = initial state)
q = zeros(narma_len, 1);

% NARMA2 recurrence:
%   q(n+1) = 0.4*q(n) + 0.4*q(n)*q(n-1) + 0.6*u(n)^3 + 0.1
for n = narma_order : num_tot_points - 1
    q(n+1) = 0.4 * q(n) + 0.4 * q(n) * q(n-1) + 0.6 * u(n)^3 + 0.1;
end
toc


%% =====================================================================
%%  OUTER LOOP — param2
%% =====================================================================
for i2 = 1:nP2
    val2 = param2_values(i2);

    % Apply param2 to the correct variable
    switch param2_name
        case 'washout1';            washout1 = val2;
        case 'num_train_points';    num_train_points = val2;
        case 'wheretostarttest'
            wheretostarttest = val2;
            washout2 = wheretostarttest - (washout1 + num_train_points);
        case 'washout2';            washout2 = val2;
        case 'num_test_points';     num_test_points = val2;
        case 'offset';              offset = val2;
        case 'lambda';              lambda = val2;
        case 'lambda_smol';         lambda_smol = val2;
        case 'movmean_window';      movmean_window = val2;
        case 'movmean_readout_test_window'; movmean_readout_test_window = val2;
        case 'Nres_init'
            Nres_init = val2;
            Nres = Nres_init * memorywindowlength;
        case 'memorywindowlength'
            memorywindowlength = val2;
            Nres = Nres_init * memorywindowlength;
        case 'memlengthsweep'
            memlengthsweep = val2;
            memory_length = round((memlengthsweep * Tpeak) / T) * T;
        otherwise
            error('Parameter2 name (%s) not recognized.', param2_name);
    end

    %% INNER LOOP — param1
    for i1 = 1:nP1
        tic
        val1 = param1_values(i1);

        % Apply param1 to the correct variable
        switch param1_name
            case 'washout1';            washout1 = val1;
            case 'num_train_points';    num_train_points = val1;
            case 'wheretostarttest'
                wheretostarttest = val1;
                washout2 = wheretostarttest - (washout1 + num_train_points);
            case 'washout2';            washout2 = val1;
            case 'num_test_points';     num_test_points = val1;
            case 'offset';              offset = val1;
            case 'lambda';              lambda = val1;
            case 'lambda_smol';         lambda_smol = val1;
            case 'movmean_window';      movmean_window = val1;
            case 'movmean_readout_test_window'; movmean_readout_test_window = val1;
            case 'Nres_init'
                Nres_init = val1;
                Nres = Nres_init * memorywindowlength;
            case 'memorywindowlength'
                memorywindowlength = val1;
                Nres = Nres_init * memorywindowlength;
            case 'memlengthsweep'
                memlengthsweep = val1;
                memory_length = round((memlengthsweep * Tpeak) / T) * T;
            otherwise
                error('Parameter1 name (%s) not recognized.', param1_name);
        end

        num_tot_points = washout1 + washout2 + num_train_points + num_test_points;

        % Molecule count mapping from the (fixed) input sequence
        N_i = N_min + (u - min(u)) / (max(u) - min(u)) * (N_max - N_min);

        %% 3. Time Vector
        t_total      = length(N_i) * T;
        t_values     = 0:dt:t_total;
        steps_per_symbol = T / dt;

        %% 4. Pre-allocate Concentration and Receptor Occupation Vectors
        c_values    = zeros(size(t_values));
        n_values    = zeros(size(t_values));
        n_values(1) = 0;

        %% 5. Compute Ligand Concentration c(t) — Deterministic Model
        for i = 1:length(N_i)
            t_symbol_start = (i - 1) * T;
            t_memory_end   = t_symbol_start + memory_length;
            indices = find(t_values > t_symbol_start & t_values <= t_memory_end);
            for j = indices'
                t = t_values(j) - t_symbol_start;
                c_values(j) = c_values(j) + (N_i(i) ./ ((4 * pi * D * t).^(3/2))) .* ...
                              exp(-distance^2 ./ (4 * D * t));
            end
        end

        %% 6. Receptor Binding ODE (Euler)
        for i = 2:length(t_values)
            c_t = c_values(i - 1);
            n_t = n_values(i - 1);
            dn_dt = k_on * (N - n_t) * c_t - k_off * n_t;
            n_values(i) = n_t + dn_dt * dt;
        end
        n_values = n_values / N;  % Normalise to [0, 1]

        if plot_in_loop
            figure('Position', [400 400 1200 400]);
            plot(t_values, n_values, 'LineWidth', 2);
            xlabel('Time (s)'); ylabel('<n(t)>');
            title('NARMA2 – Bound Receptors (Deterministic)'); grid on;

            figure('Position', [400 400 1200 400]);
            plot(t_values, c_values, 'LineWidth', 2);
            xlabel('Time (s)'); ylabel('c(t)');
            title('NARMA2 – Ligand Concentration'); grid on;
        end

        %% 7. Training Phase — Deterministic
        sampled_matrix             = zeros(Nres, num_train_points);
        sample_indices_to_be_deleted = [];

        for i = washout1 + 1 : washout1 + num_train_points
            t_symbol_start = (i - 1) * T + offset;
            start_idx      = find(t_values >= t_symbol_start, 1);
            sample_indices = start_idx - (memorywindowlength - 1) * steps_per_symbol + ...
                             (0:Nres - 1) * (steps_per_symbol / Nres) * memorywindowlength;
            sample_indices(sample_indices > length(t_values)) = [];

            if length(sample_indices) >= Nres
                sampled_matrix(:, i - washout1) = n_values(sample_indices(1:Nres));
            else
                sample_indices_to_be_deleted = [sample_indices_to_be_deleted, i - washout1]; %#ok<AGROW>
            end
        end
        sampled_matrix(:, sample_indices_to_be_deleted) = [];

        % Bias term
        sampled_matrix = [sampled_matrix; ones(1, size(sampled_matrix, 2))];

        % Target: q(n+1) aligned to training window
        % (same indexing as NARMA10 sweep — +2 shift accounts for NARMA warmup)
        q_train = q(washout1 + 2 : washout1 + num_train_points + 1)';
        q_train = q_train(:);

        sampled_matrix = sampled_matrix';
        W_out = pinv(sampled_matrix' * sampled_matrix + lambda * eye(size(sampled_matrix, 2))) ...
                * (sampled_matrix' * q_train);

        q_train_hat = sampled_matrix * W_out;

        if plot_in_loop
            figure('Position', [400 400 1200 400]);
            plot(q_train(1:end-2),     'b-o', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(q_train_hat(1:end-2), 'r-',  'LineWidth', 2);
            xlabel('Index'); ylabel('Target / Estimate');
            title('NARMA2 – Training: Deterministic'); legend('Target','Estimate'); grid on;
        end

        %% 8. NRMSE — Training (Deterministic)
        rmse_train  = sqrt(mean((q_train_hat(1:end-2) - q_train(1:end-2)).^2));
        nrmse_train = rmse_train / std(q_train(1:end-2));

        %% 9. Testing Phase — Deterministic
        sampled_matrix_test          = zeros(Nres, num_test_points);
        sample_indices_to_be_deleted = [];

        for i = washout1 + num_train_points + washout2 + 1 : num_tot_points
            t_symbol_start = (i - 1) * T + offset;
            start_idx      = find(t_values >= t_symbol_start, 1);
            sample_indices = start_idx - (memorywindowlength - 1) * steps_per_symbol + ...
                             (0:Nres - 1) * (steps_per_symbol / Nres) * memorywindowlength;
            sample_indices(sample_indices > length(t_values)) = [];

            if length(sample_indices) >= Nres
                sampled_matrix_test(:, i - (washout1 + num_train_points + washout2)) = ...
                    n_values(sample_indices(1:Nres));
            else
                sample_indices_to_be_deleted = [sample_indices_to_be_deleted, ...
                    i - (washout1 + num_train_points + washout2)]; %#ok<AGROW>
            end
        end
        sampled_matrix_test(:, sample_indices_to_be_deleted) = [];
        sampled_matrix_test = [sampled_matrix_test; ones(1, size(sampled_matrix_test, 2))];
        sampled_matrix_test = sampled_matrix_test';

        % Target for test window
        q_test = q(washout1 + num_train_points + washout2 + 2 : num_tot_points + 1)';
        q_test = q_test(:);

        q_test_hat = sampled_matrix_test * W_out;

        %% 10. NRMSE — Testing (Deterministic)
        rmse_test  = sqrt(mean((q_test_hat(1:end-2) - q_test(1:end-2)).^2));
        nrmse_test = rmse_test / std(q_test(1:end-2));

        if plot_in_loop
            figure('Position', [400 400 1200 400]);
            plot(q_test(1:end-2),     'b-o', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(q_test_hat(1:end-2), 'r-',  'LineWidth', 2);
            xlabel('Index'); ylabel('Target / Estimate');
            title('NARMA2 – Testing: Deterministic'); legend('Target','Estimate'); grid on;
        end

        %% Save deterministic results
        nrmse_train_NUM_array(i1, i2) = nrmse_train;
        nrmse_test_NUM_array(i1, i2)  = nrmse_test;

        fprintf('%s=%g, %s=%g -> NRMSE-NUM(Tr)=%.3g, NRMSE-NUM(Te)=%.3g\n', ...
            param1_name, val1, param2_name, val2, nrmse_train, nrmse_test);

        if smold == 0
            toc
            continue
        end

        %% ============================================================
        %%  SMOLDYN SECTION
        %% ============================================================

        numIterations = 1;

        KON               = k_on  * 1e18;
        KOFF              = k_off;
        DIFFMESSENGER     = D     * 1e12;
        DIFFRECEPTOR      = 0;
        DIFFRECEPTORACT   = 0;
        START_TIME        = 0;
        BIT_INTERVAL      = T;
        STOP_TIME         = START_TIME + length(N_i) * BIT_INTERVAL;
        TIME_STEP         = 0.01;
        SAMPLING_PERIOD   = 1;
        numRECEPTOR       = N;
        TXRXDISTANCE      = distance * 1e6;
        BOUNDARYLENGTH    = 25;
        TXPOSITION        = -10;
        RXRADIUS          = 3;
        RXRECEPTIONSPACETHICKNESS = 0.2;

        % Build Smoldyn config (generateSimConfig is unchanged)
        input_time_points = START_TIME : BIT_INTERVAL : STOP_TIME - eps * 1e10;
        % generateSimConfig(N_i, input_time_points)  % Uncomment to regenerate

        % Directory for this run (or use existing data dir)
        dirname = sprintf('NARMA2_sample_%s/', datestr(now, 'mm-dd-yyyy--HH-MM-SS'));
        % mkdir(dirname)
        % copyfile('MCpointTXsphRX_generated.txt', dirname);

        variableNameArray = [...
            "KON","KOFF","DIFFMESSENGER","DIFFRECEPTOR","DIFFRECEPTORACT",...
            "START_TIME","STOP_TIME","TIME_STEP","BIT_INTERVAL","numRECEPTOR",...
            "TXRXDISTANCE","BOUNDARYLENGTH","TXPOSITION","RXRADIUS",...
            "RXRECEPTIONSPACETHICKNESS","SAMPLING_PERIOD"];

        variableNo     = 2;
        variableValues = KOFF;

        % -- Run Smoldyn (skipped when smold_only_data = 1) --
        for j = 1:numel(variableValues)
            assignin('base', variableNameArray(variableNo), variableValues(j));
            for ii = 1:numIterations
                command = sprintf(['smoldyn %sMCpointTXsphRX_generated.txt -wt -s ', ...
                    '--define index=%d --define variableNo=%d --define variableValue=%d ',...
                    '--define KON=%d --define KOFF=%d --define DIFFMESSENGER=%d ',...
                    '--define DIFFRECEPTOR=%d --define DIFFRECEPTORACT=%d ',...
                    '--define START_TIME=%d --define STOP_TIME=%d --define TIME_STEP=%d ',...
                    '--define BIT_INTERVAL=%d --define numRECEPTOR=%d ',...
                    '--define TXRXDISTANCE=%d --define BOUNDARYLENGTH=%d ',...
                    '--define TXPOSITION=%d --define RXRADIUS=%d ',...
                    '--define RXRECEPTIONSPACETHICKNESS=%d --define SAMPLING_PERIOD=%d '], ...
                    dirname, ii, variableNo, variableValues(j), ...
                    KON, KOFF, DIFFMESSENGER, DIFFRECEPTOR, DIFFRECEPTORACT, ...
                    START_TIME, STOP_TIME, TIME_STEP, BIT_INTERVAL, numRECEPTOR, ...
                    TXRXDISTANCE, BOUNDARYLENGTH, TXPOSITION, RXRADIUS, ...
                    RXRECEPTIONSPACETHICKNESS, SAMPLING_PERIOD);

                if smold_only_data ~= 1
                    system(command);
                end
            end
        end

        % -- Read Smoldyn output files --
        if smold_only_data == 1
            dirname = sprintf('%s/', datadirname);
        end

        received_signals = cell(numel(variableValues), numIterations);
        times1           = cell(numel(variableValues), numIterations);
        ReceptorActives  = cell(numel(variableValues), numIterations);
        times2           = cell(numel(variableValues), numIterations);

        for j = 1:numel(variableValues)
            assignin('base', variableNameArray(variableNo), variableValues(j));
            for k = 1:numIterations
                receivedSignalFilename = sprintf('%sreceivedsignal_varNo_%d_varValue_%d_iter_%d_NARMA2.txt', ...
                    dirname, variableNo, variableValues(j), k);
                data1 = importdata(receivedSignalFilename, ' ');
                received_signals{j, k} = data1(:, 3);
                times1{j, k}           = data1(:, 1);

                allMoleculesFilename = sprintf('%sallmolecules_varNo_%d_varValue_%d_iter_%d.txt', ...
                    dirname, variableNo, variableValues(j), k);
                data2_temp           = importdata(allMoleculesFilename, ' ', 1);
                data2                = data2_temp.data;
                times2{j, k}         = data2(:, 1);
                ReceptorActives{j, k} = data2(:, end);
            end
        end

        % Extract and optionally filter the active receptor signal
        received_signal_fin = received_signals{1, 1};
        Active_rec_fin      = ReceptorActives{1, 1} / N;

        if movmean_on == 1
            Active_rec_fin = movmean(Active_rec_fin, [movmean_window 0]); % causal average
        end

        if plot_in_loop
            figure('Position', [400 400 1200 400]);
            plot(t_values, n_values, 'LineWidth', 2); hold on;
            plot(linspace(START_TIME, STOP_TIME, length(Active_rec_fin)), Active_rec_fin, '-');
            grid on; xlabel('Time (s)'); ylabel('<n(t)>');
            title('NARMA2 – Bound Receptors: Deterministic vs. Smoldyn');
            legend('Deterministic','Smoldyn');
        end

        %% 11. Training Phase — Smoldyn
        sampled_matrix_smol          = zeros(Nres, num_train_points);
        time_index_smol              = START_TIME : TIME_STEP * SAMPLING_PERIOD : STOP_TIME;
        steps_per_symbol_smol        = T / TIME_STEP;
        sample_indices_to_be_deleted = [];

        for i = washout1 + 1 : washout1 + num_train_points
            t_symbol_start_smol = (i - 1) * T + offset;
            start_idx           = find(time_index_smol >= t_symbol_start_smol, 1);
            sample_indices_smol = start_idx - (memorywindowlength - 1) * steps_per_symbol_smol + ...
                                  (0:Nres - 1) * (steps_per_symbol_smol / Nres) * memorywindowlength;
            sample_indices_smol(sample_indices_smol > length(time_index_smol)) = [];

            if length(sample_indices_smol) >= Nres
                sampled_matrix_smol(:, i - washout1) = Active_rec_fin(sample_indices_smol(1:Nres));
            else
                sample_indices_to_be_deleted = [sample_indices_to_be_deleted, i - washout1]; %#ok<AGROW>
            end
        end
        sampled_matrix_smol(:, sample_indices_to_be_deleted) = [];
        sampled_matrix_smol = [sampled_matrix_smol; ones(1, size(sampled_matrix_smol, 2))];
        sampled_matrix_smol = sampled_matrix_smol';

        W_out_smol = pinv(sampled_matrix_smol' * sampled_matrix_smol + ...
                          lambda_smol * eye(size(sampled_matrix_smol, 2))) ...
                     * (sampled_matrix_smol' * q_train);

        q_train_hat_smol = sampled_matrix_smol * W_out_smol;

        if plot_in_loop
            figure('Position', [400 400 1200 400]);
            plot(q_train(1:end-2),          'b-o', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(q_train_hat_smol(1:end-2), 'r-',  'LineWidth', 2);
            xlabel('Index'); ylabel('Target / Estimate');
            title('NARMA2 – Training: Smoldyn'); legend('Target','Estimate'); grid on;
        end

        %% 12. NRMSE — Training (Smoldyn)
        q_train          = q_train(:);
        q_train_hat_smol = q_train_hat_smol(:);

        rmse_train_smol  = sqrt(mean((q_train_hat_smol(1:end-2) - q_train(1:end-2)).^2));
        std_q_train_smol = std(q_train(1:end-2));
        nrmse_train_smol = rmse_train_smol / std_q_train_smol;

        %% 13. Testing Phase — Smoldyn
        sampled_matrix_test_smol     = zeros(Nres, num_test_points);
        sample_indices_to_be_deleted = [];

        for i = washout1 + num_train_points + washout2 + 1 : num_tot_points
            t_symbol_start_smol = (i - 1) * T + offset;
            start_idx           = find(time_index_smol >= t_symbol_start_smol, 1);
            sample_indices_smol = start_idx - (memorywindowlength - 1) * steps_per_symbol_smol + ...
                                  (0:Nres - 1) * (steps_per_symbol_smol / Nres) * memorywindowlength;
            sample_indices_smol(sample_indices_smol > length(time_index_smol)) = [];

            if length(sample_indices_smol) >= Nres
                sampled_matrix_test_smol(:, i - (washout1 + num_train_points + washout2)) = ...
                    Active_rec_fin(sample_indices_smol(1:Nres));
            else
                sample_indices_to_be_deleted = [sample_indices_to_be_deleted, ...
                    i - (washout1 + num_train_points + washout2)]; %#ok<AGROW>
            end
        end
        sampled_matrix_test_smol(:, sample_indices_to_be_deleted) = [];
        sampled_matrix_test_smol = [sampled_matrix_test_smol; ones(1, size(sampled_matrix_test_smol, 2))];
        sampled_matrix_test_smol = sampled_matrix_test_smol';

        q_test_hat_smol = sampled_matrix_test_smol * W_out_smol;

        %% 14. NRMSE — Testing (Smoldyn)
        if movmean_readout_test == 1
            q_test_hat_smol = movmean(q_test_hat_smol, movmean_readout_test_window);
        end

        rmse_test_smol  = sqrt(mean((q_test_hat_smol(1:end-2) - q_test(1:end-2)).^2));
        std_q_test_smol = std(q_test(1:end-2));
        nrmse_test_smol = rmse_test_smol / std_q_test_smol;

        % Final test plot (always shown, mirrors NARMA10 sweep)
        figure('Position', [400 400 1200 400]);
        plot(q_test(1:end-2),          'b-o', 'LineWidth', 2, 'MarkerSize', 6);
        hold on;
        plot(q_test_hat_smol(1:end-2), 'r-',  'LineWidth', 2);
        xlabel('Index'); ylabel('Target / Estimate');
        title(sprintf('NARMA2 – Testing: Smoldyn  |  NRMSE = %.4f', nrmse_test_smol));
        legend('Target','Smoldyn Estimate'); grid on;

        %% Save Smoldyn results
        nrmse_train_SMOL_array(i1, i2) = nrmse_train_smol;
        nrmse_test_SMOL_array(i1, i2)  = nrmse_test_smol;

        fprintf('%s=%g, %s=%g -> NRMSE-SMOL(Tr)=%.3g, NRMSE-SMOL(Te)=%.3g\n', ...
            param1_name, val1, param2_name, val2, nrmse_train_smol, nrmse_test_smol);
        toc

    end % END inner loop (param1)
end % END outer loop (param2)


%% =====================================================================
%%  2-D HEATMAPS
%% =====================================================================
[X, Y]      = meshgrid(param1_values, param2_values);
figureOpts  = {'LineColor','none','LevelStepMode','auto'};

% Deterministic NRMSE
figure('Name','NARMA2 – NRMSE Heatmap (Deterministic)');
contourf(X, Y, nrmse_test_NUM_array', 20, figureOpts{:});
shading interp; colorbar; colormap(flipud(parula));
if useLogX, set(gca,'XScale','log'); end
if useLogY, set(gca,'YScale','log'); end
xlabel(param1_name); ylabel(param2_name);
title('NARMA2 – Test NRMSE Heatmap (Deterministic Model)');

% Smoldyn NRMSE
figure('Name','NARMA2 – NRMSE Heatmap (Stochastic / Smoldyn)');
contourf(X, Y, nrmse_test_SMOL_array', 20, figureOpts{:});
shading interp; colorbar; colormap(flipud(parula));
if useLogX, set(gca,'XScale','log'); end
if useLogY, set(gca,'YScale','log'); end
xlabel(param1_name); ylabel(param2_name);
title('NARMA2 – Test NRMSE Heatmap (Stochastic / Smoldyn)');


%% Appendix
% (Mac users) Add Smoldyn to PATH if MATLAB cannot find it:
%   setenv('PATH', getenv('PATH') + ":/usr/local/bin")
%   system('echo $PATH')
