%%
% Code to generate Smoldyn data given the MC system settings and input
% sequence generated for NARMA2 series.
%
% NARMA2 recurrence (order-2 NARMA):
%   q(n+1) = 0.4*q(n) + 0.4*q(n)*q(n-1) + 0.6*u(n)^3 + 0.1
%
%%

%% Initialize the MATLAB environment
close all   % Close all open figure windows
clear       % Clear all variables from the workspace

movmean_on = 0; % = 0 by default
movmean_window = 500; % only meaningful if movmean_on = 1

rng(0); % IMPORTANT TO FIX THE RANDOM SEQUENCE (both for numerical and Smoldyn)

%% 1. Parameter Setup
% Receptor and binding parameters
N = 500;                % Total number of receptors at the receiver
k_on = 1e-18;          % Binding rate constant (in m^3/(mol*s), i.e. 1/(M*s) scaled)
k_off = 1;             % Unbinding rate constant (in 1/s)
KD = k_off / k_on;    % Dissociation constant

% Time parameters
T = 1;                   % Symbol duration (seconds)

washout1 = 100;           % Number of symbols to wash out before training
num_train_points = 2000;  % Number of symbols used for training
wheretostarttest = 2500;
washout2 = wheretostarttest - (washout1 + num_train_points); % Symbols to wash out before testing
num_test_points = 2500;   % Number of symbols used for testing

num_tot_points = washout1 + washout2 + num_train_points + num_test_points; % Total number of symbols

% Communication channel parameters
distance = 10e-6;        % Distance between transmitter and receiver (10 microns)
D = 1e-11;               % Diffusion coefficient (in m^2/s)

% Input normalization parameters
N_min = 100;             % Minimum value for input normalization
N_max = 3000;            % Maximum value for input normalization

% Step size control
dt = 0.001;              % Time step for deterministic/numerical analysis (fixed)
memlengthsweep = 100;
Tpeak = distance^2 / (6 * D);
memory_length = round((memlengthsweep * Tpeak) / T) * T; % memory of MC channel
offset = 0;


%% 2. Generate the NARMA2 Time Series (Only ONCE)

% NARMA2 requires a warmup of only 2 samples (lag-2 system).
% We generate (num_tot_points + 2 + 1) values:
%   +2  for warmup (indices 1..2 are initialised to zero)
%   +1  so we have next-step target pairs throughout
narma_order = 2;
narma_len   = num_tot_points + narma_order + 1;

% Random input in [0, 0.5]
u = 0.5 * rand(narma_len, 1);

% Pre-allocate the NARMA2 output (initialised to zero)
q = zeros(narma_len, 1);

% Compute q(n) for the NARMA2 system.
% Standard NARMA2 formula:
%   q(n+1) = 0.4*q(n) + 0.4*q(n)*q(n-1) + 0.6*u(n)^3 + 0.1
% The loop starts at n = narma_order so that q(n-1) is always defined.
for n = narma_order : num_tot_points - 1
    q(n+1) = 0.4 * q(n) + 0.4 * q(n) * q(n-1) + 0.6 * u(n)^3 + 0.1;
end

% Normalize input u(n) to obtain molecule counts N_i
N_i = N_min + (u - min(u)) / (max(u) - min(u)) * (N_max - N_min);


%% 3. Generate Time Vector
t_total = length(N_i) * T;      % Total simulation time
t_values = 0:dt:t_total;        % Time vector from 0 to t_total with step dt

% Number of time steps per symbol duration
steps_per_symbol = T / dt;


%% 4. Pre-allocate Concentration and Receptor Occupation Vectors
c_values = zeros(size(t_values));  % Ligand concentration at the receiver
n_values = zeros(size(t_values));  % Bound-receptor count (normalised later)

% Initial condition: no bound receptors at t = 0
n_values(1) = 0;


%% 5. Compute Ligand Concentration c(t) at the Receiver
tic
% Superpose the free-space 3-D diffusion Green's function for each symbol.
for i = 1:length(N_i)
    t_symbol_start = (i - 1) * T;             % Start time of current symbol
    t_memory_end   = t_symbol_start + memory_length; % End of memory window

    % Indices within the memory window
    indices = find(t_values > t_symbol_start & t_values <= t_memory_end);

    for j = indices'
        t = t_values(j) - t_symbol_start;    % Time elapsed since transmission
        c_values(j) = c_values(j) + (N_i(i) ./ ((4 * pi * D * t).^(3/2))) .* ...
                      exp(-distance^2 ./ (4 * D * t));
    end
end
toc


%% 6. Simulate Receptor Binding Process (Euler's Method)
for i = 2:length(t_values)
    c_t   = c_values(i - 1);  % Concentration at previous time step
    n_t   = n_values(i - 1);  % Bound receptors at previous time step

    % Mean-field receptor binding ODE
    dn_dt = k_on * (N - n_t) * c_t - k_off * n_t;

    n_values(i) = n_t + dn_dt * dt;
end

% Normalise to average occupation ratio in [0, 1]
n_values = n_values / N;

% --- Plots for deterministic model ---
figure('Position', [400 400 1200 400]);
plot(t_values, n_values, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('<n(t)> - Average Fraction of Bound Receptors');
title('NARMA2 – Average Fraction of Bound Receptors vs. Time (Deterministic)');
grid on;

figure('Position', [400 400 1200 400]);
plot(t_values, c_values, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('c(t) - Ligand Concentration');
title('NARMA2 – Ligand Concentration vs. Time at the Receiver');
grid on;


%% =====================================================================
%%  SMOLDYN SECTION
%% =====================================================================

tic

%% Smoldyn Simulation Parameters
numIterations = 1;  % Number of Monte Carlo iterations

% Reaction rates / diffusion (Smoldyn unit conventions)
KON                    = k_on * 1e18;  % On rate  [um^3 / (mol*s)] -> Smoldyn units
KOFF                   = k_off;        % Off rate [1/s]
DIFFMESSENGER          = D * 1e12;     % Diffusion coefficient [um^2/s]
DIFFRECEPTOR           = 0;
DIFFRECEPTORACT        = 0;

% Simulation time parameters
START_TIME             = 0;
BIT_INTERVAL           = T;
STOP_TIME              = START_TIME + length(N_i) * BIT_INTERVAL;
disp(['Smoldyn Simulation Stop Time: ', num2str(STOP_TIME)]);
TIME_STEP              = 0.01;
SAMPLING_PERIOD        = 1;   % In time steps

% Receiver
numRECEPTOR            = N;

% Geometry
TXRXDISTANCE           = distance * 1e6;  % [um]
BOUNDARYLENGTH         = 25;
TXPOSITION             = -10;
RXRADIUS               = 3;
RXRECEPTIONSPACETHICKNESS = 0.2;


%% Create Smoldyn Configuration File
% generateSimConfig() is unchanged — it only depends on N_i and the
% time points, not on the specific NARMA order.
input_time_points = START_TIME : BIT_INTERVAL : STOP_TIME - eps * 1e10;
generateSimConfig(N_i, input_time_points)


%% Create Time-Stamped Output Directory
dirname = sprintf('NARMA2_sample_%s/', datestr(now, 'mm-dd-yyyy--HH-MM-SS'));
mkdir(dirname)

status = copyfile('MCpointTXsphRX_generated.txt', dirname);
if ~status
    disp('Smoldyn execution file could not be copied to new directory.')
    return
end


%% Define Variables for Parameter Sweep
variableNameArray = [...
    "KON",...                            % 1
    "KOFF",...                           % 2
    "DIFFMESSENGER",...                  % 3
    "DIFFRECEPTOR",...                   % 4
    "DIFFRECEPTORACT",...                % 5
    "START_TIME",...                     % 6
    "STOP_TIME",...                      % 7
    "TIME_STEP",...                      % 8
    "BIT_INTERVAL",...                   % 9
    "numRECEPTOR",...                    % 10
    "TXRXDISTANCE",...                   % 11
    "BOUNDARYLENGTH",...                 % 12
    "TXPOSITION",...                     % 13
    "RXRADIUS",...                       % 14
    "RXRECEPTIONSPACETHICKNESS",...      % 15
    "SAMPLING_PERIOD"...                 % 16
    ];

variableNo     = 2;       % Sweep KOFF (index 2)
variableValues = KOFF;    % Single value — no sweep, just run at default


%% Run Smoldyn Simulations
for j = 1:numel(variableValues)
    assignin('base', variableNameArray(variableNo), variableValues(j));

    for i = 1:numIterations
        command = sprintf(['smoldyn %sMCpointTXsphRX_generated.txt -wt -s ', ...
            '--define index=%d ',...
            '--define variableNo=%d ',...
            '--define variableValue=%d ',...
            '--define KON=%d ',...
            '--define KOFF=%d ',...
            '--define DIFFMESSENGER=%d ',...
            '--define DIFFRECEPTOR=%d ',...
            '--define DIFFRECEPTORACT=%d ',...
            '--define START_TIME=%d ',...
            '--define STOP_TIME=%d ',...
            '--define TIME_STEP=%d ',...
            '--define BIT_INTERVAL=%d ',...
            '--define numRECEPTOR=%d ',...
            '--define TXRXDISTANCE=%d ',...
            '--define BOUNDARYLENGTH=%d ',...
            '--define TXPOSITION=%d ',...
            '--define RXRADIUS=%d ',...
            '--define RXRECEPTIONSPACETHICKNESS=%d ',...
            '--define SAMPLING_PERIOD=%d '], ...
            dirname, i, variableNo, variableValues(j), ...
            KON, KOFF, DIFFMESSENGER, DIFFRECEPTOR, DIFFRECEPTORACT, ...
            START_TIME, STOP_TIME, TIME_STEP, BIT_INTERVAL, numRECEPTOR, ...
            TXRXDISTANCE, BOUNDARYLENGTH, TXPOSITION, RXRADIUS, ...
            RXRECEPTIONSPACETHICKNESS, SAMPLING_PERIOD);

        system(command);
    end
end


%% Read Simulation Results
received_signals = cell(numel(variableValues), numIterations);
times1           = cell(numel(variableValues), numIterations);
ReceptorActives  = cell(numel(variableValues), numIterations);
times2           = cell(numel(variableValues), numIterations);
bitsequences     = cell(numel(variableValues), numIterations);

for j = 1:numel(variableValues)
    assignin('base', variableNameArray(variableNo), variableValues(j));

    for k = 1:numIterations
        fprintf('varName = %s -- varValue = %g -- numIter = %d\n', ...
            variableNameArray(variableNo), variableValues(j), k);

        receivedSignalFilename = sprintf('%sreceivedsignal_varNo_%d_varValue_%d_iter_%d.txt', ...
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


%% Post-process and Plot Smoldyn Results
received_signal_fin = received_signals{1, 1};
Active_rec_fin      = ReceptorActives{1, 1} / N;

if movmean_on == 1
    % Causal moving average (past values only) — same convention as stochastic sweep
    Active_rec_fin = movmean(Active_rec_fin, [movmean_window 0]);
end

figure('Position', [400 400 1200 400]);
plot(t_values, n_values, 'LineWidth', 2);
hold on
plot(linspace(START_TIME, STOP_TIME, length(Active_rec_fin)), Active_rec_fin, '-');
grid on;
xlabel('Time (s)');
ylabel('<n(t)> - Average Fraction of Bound Receptors');
title('NARMA2 – Bound Receptors: Deterministic vs. Smoldyn');
legend('Deterministic', 'Smoldyn');

figure('Position', [400 400 1200 400]);
plot(linspace(START_TIME, STOP_TIME, length(received_signal_fin)), received_signal_fin, '-');
grid on;
xlabel('Time (s)');
ylabel('Ligands in reception space (Smoldyn)');
title('NARMA2 – Received Signal (Smoldyn)');

toc


%% Appendix
% (Mac users) Add Smoldyn to PATH if MATLAB cannot find it:
%   setenv('PATH', getenv('PATH') + ":/usr/local/bin")
%   system('echo $PATH')
