%% MonteCarloConfirmation.m
% Monte Carlo validation for point-source / point-receiver approximation
%
% Runs:
%   1) Deterministic analytical point model
%   2) Smoldyn point-TX / finite spherical RX
%   3) Smoldyn finite-TX / finite spherical RX
%
% Input: 10 pulses spaced by T = 1 s

%% 0a. Initialize MATLAB
close all
clear
clc

%% 0b. Control flags

smold = 1;                 % 1 = run/read Smoldyn
smold_only_data = 0;       % 1 = only read existing data
datadirname = '';          % only used if smold_only_data = 1

smoldyn_visualization_on = 0;
% 0 = no graphics
% 1 = graphics on

%% 1. Parameter Setup

% Receptor and binding parameters
N = 500;                % Total number of receptors at receiver
k_on = 1.e-18;            % Binding rate constant, m^3/s
k_off = 1;               % Unbinding rate constant, 1/s
KD = k_off / k_on;

% Time parameters
T = 1;                   % Pulse spacing / symbol duration
numPulses = 20;
START_TIME = 0;
STOP_TIME = 25;          % Total simulation runtime, seconds

% Communication channel parameters
distance = 10e-6;        % 10 microns
D = 1e-11;               % m^2/s

% Input molecule count range
N_min = 100;
N_max = 3000;

% 10-pulse input
N_i = round(0.5 * (N_min + N_max)) * ones(1, numPulses);
input_time_points = START_TIME:T:(numPulses-1)*T;

% Smoldyn simulation parameters
numIterations = 20;      % Monte Carlo trials per stochastic configuration
TIME_STEP = 0.001;       % 1 ms
SAMPLING_PERIOD = 5;     % Output every 5 time steps = 5 ms

% Smoldyn reaction/diffusion parameters
KON = k_on * 1e18;
KOFF = k_off;
DIFFMESSENGER = D * 1e12;
DIFFRECEPTOR = 0;
DIFFRECEPTORACT = 0;

BIT_INTERVAL = T;
numRECEPTOR = N;

% Spatial parameters (Smoldyn units are micrometers)
TXRXDISTANCE = distance * 1e6;
BOUNDARYLENGTH = 30;
TXPOSITION = -10;
RXRADIUS = 3;                       
RXRECEPTIONSPACETHICKNESS = 0.2;

% Finite transmitter radii to test
finiteTxRadii = [0.5];    

% Smoldyn executable
smoldynCmd = 'smoldyn';

if smoldyn_visualization_on == 1
    smoldyn_flag = '-w -s';
else
    smoldyn_flag = '-wt -s';
end

%% 2. Deterministic analytical point model

dt_det = TIME_STEP * SAMPLING_PERIOD;

[t_det, b_det] = deterministicPointModel( ...
    N_i, input_time_points, STOP_TIME, dt_det, distance, D, k_on, k_off);

figure;
plot(t_det, b_det, 'k-', 'LineWidth', 2)
grid on
xlabel('Time (s)')
ylabel('b(t)')
title('Deterministic analytical point model')

if smold == 0
    return
end

%% 3. Define Smoldyn variable list

variableNo = 2;
variableValues = KOFF;

%% 4. Run/read case 1: Smoldyn point-TX / finite spherical RX

case_point.dirname = sprintf('sample_pointTX_finiteRX_%s/', datestr(now,'mm-dd-yyyy--HH-MM-SS'));
case_point.txMode = "point";
case_point.r_tx = 0;

[t_point, B_point] = runOneCase( ...
    case_point, ...
    N_i, input_time_points, ...
    smoldyn_visualization_on, smoldyn_flag, smoldynCmd, ...
    smold_only_data, datadirname, ...
    variableNo, variableValues, ...
    KON, KOFF, DIFFMESSENGER, DIFFRECEPTOR, DIFFRECEPTORACT, ...
    START_TIME, STOP_TIME, TIME_STEP, BIT_INTERVAL, numRECEPTOR, ...
    TXRXDISTANCE, BOUNDARYLENGTH, TXPOSITION, RXRADIUS, ...
    RXRECEPTIONSPACETHICKNESS, SAMPLING_PERIOD, N, numIterations);

b_point_mean = mean(B_point, 2, 'omitnan');
b_point_std = std(B_point, 0, 2, 'omitnan');

%% 5. Run/read case 2: Smoldyn finite-TX / finite spherical RX

finiteResults = struct([]);

for rr = 1:length(finiteTxRadii)

    r_tx = finiteTxRadii(rr);

    case_finite.dirname = sprintf('sample_finiteTX_r%.2fum_finiteRX_%s/', ...
        r_tx, datestr(now,'mm-dd-yyyy--HH-MM-SS'));
    case_finite.txMode = "finite";
    case_finite.r_tx = r_tx;

    [t_fin, B_fin] = runOneCase( ...
        case_finite, ...
        N_i, input_time_points, ...
        smoldyn_visualization_on, smoldyn_flag, smoldynCmd, ...
        smold_only_data, datadirname, ...
        variableNo, variableValues, ...
        KON, KOFF, DIFFMESSENGER, DIFFRECEPTOR, DIFFRECEPTORACT, ...
        START_TIME, STOP_TIME, TIME_STEP, BIT_INTERVAL, numRECEPTOR, ...
        TXRXDISTANCE, BOUNDARYLENGTH, TXPOSITION, RXRADIUS, ...
        RXRECEPTIONSPACETHICKNESS, SAMPLING_PERIOD, N, numIterations);

    finiteResults(rr).r_tx = r_tx;
    finiteResults(rr).t = t_fin;
    finiteResults(rr).B = B_fin;
    finiteResults(rr).mean = mean(B_fin, 2, 'omitnan');
    finiteResults(rr).std = std(B_fin, 0, 2, 'omitnan');
end

%% 6. Metrics

fprintf('\n================ Metrics vs deterministic analytical point model ================\n')

b_det_point = interp1(t_det, b_det, t_point, 'linear', 'extrap');
metrics_point = computeTraceMetrics(b_det_point, b_point_mean, t_point);

fprintf('Point-TX / finite-RX: corr = %.4f, NRMSE = %.4f, peak-time error = %.4f s\n', ...
    metrics_point.corr, metrics_point.nrmse, metrics_point.peakTimeError)

for rr = 1:length(finiteResults)

    b_det_fin = interp1(t_det, b_det, finiteResults(rr).t, 'linear', 'extrap');
    metrics_fin = computeTraceMetrics(b_det_fin, finiteResults(rr).mean, finiteResults(rr).t);
    finiteResults(rr).metrics = metrics_fin;

    fprintf('Finite-TX r_tx = %.2f um / finite-RX: corr = %.4f, NRMSE = %.4f, peak-time error = %.4f s\n', ...
        finiteResults(rr).r_tx, metrics_fin.corr, metrics_fin.nrmse, metrics_fin.peakTimeError)
end

%% 7. Final plot

figure('Color','w')
hold on
grid on
box on

plot(t_det, b_det, 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Analytical point model')

plotWithShade(t_point, b_point_mean, b_point_std, ...
    'Smoldyn point-TX / finite-RX mean')

for rr = 1:length(finiteResults)
    plotWithShade(finiteResults(rr).t, finiteResults(rr).mean, finiteResults(rr).std, ...
        sprintf('Smoldyn finite-TX r_{tx}=%.1f um / finite-RX mean', finiteResults(rr).r_tx))
end

xlabel('Time (s)')
ylabel('Receptor occupancy ratio b(t)')
title('Point approximation validation with 10-pulse input')
legend('Location','best')

savefig('pointApproxValidation.fig')
saveas(gcf, 'pointApproxValidation.png')

save('pointApproxValidation_results.mat')

fprintf('\nDone. Saved pointApproxValidation.fig/png and pointApproxValidation_results.mat\n')