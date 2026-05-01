clear; clc; close all;

%% Paper-consistent parameters
filename = 'receivedsignal_varNo_2_varValue_1.000000e-01_iter_1.txt';
T = 1.0;          % symbol duration [s]                     %NARMA10: receivedsignal_varNo_2_varValue_1_iter_1
M = 100;          % number of virtual nodes                 %MG: receivedsignal_varNo_2_varValue_1.000000e-01_iter_1
dt_expected = T/M;   % 0.01 s                               
                                                            
sigCol = 3;       % third column = received signal

%% Load data
data = readmatrix(filename);
t = data(:,1);
y = data(:,sigCol);

good = ~(isnan(t) | isnan(y));
t = t(good);
y = y(good);

%% Check sampling consistency
dt = median(diff(t));
fprintf('Median dt from file = %.6f s\n', dt);
fprintf('Expected dt = %.6f s\n', dt_expected);

if abs(dt - dt_expected) > 1e-6
    warning('File sampling interval does not match T/M exactly.');
end

%% Keep only full symbol intervals
samplesPerSymbol = M;   % since dt = T/M
numFullSymbols = floor(length(y) / samplesPerSymbol);

y = y(1:numFullSymbols * samplesPerSymbol);
t = t(1:numFullSymbols * samplesPerSymbol);

%% Build state matrix X
% each row = one symbol, each column = one virtual node
X = reshape(y, samplesPerSymbol, numFullSymbols).';
% size(X) = [numSymbols x 100]

%% Optional: discard initial transient symbols
numDiscard = 5;
if size(X,1) > numDiscard
    X = X(numDiscard+1:end,:);
end

%% Correlation / covariance
Xc = X - mean(X,1);
C = (Xc' * Xc) / size(Xc,1);
R = corrcoef(X);

%% Summary metrics
avgOffDiagAbs = (sum(abs(R),'all') - sum(abs(diag(R)))) / (M^2 - M);

corrVsSep = zeros(M-1,1);
for k = 1:M-1
    corrVsSep(k) = mean(diag(R,k));
end

eigVals = sort(eig(C), 'descend');
participationRatio = (sum(eigVals)^2) / sum(eigVals.^2);

p = eigVals / sum(eigVals);
effectiveRank = exp(-sum(p .* log(p + eps)));

fprintf('Usable symbols: %d\n', size(X,1));
fprintf('Average |off-diagonal correlation| = %.6f\n', avgOffDiagAbs);
fprintf('Participation ratio = %.4f\n', participationRatio);
fprintf('Effective rank = %.4f\n', effectiveRank);

%% Plots
figure;
imagesc(R);
axis image;
colorbar;
xlabel('Node index');
ylabel('Node index');
title('Node correlation matrix');

figure;
plot(1:M-1, corrVsSep, 'o-', 'LineWidth', 1.2);
grid on;
xlabel('\Delta node index');
ylabel('Mean correlation');
title('Correlation vs node separation');

figure;
plot(eigVals, 'o-', 'LineWidth', 1.2);
grid on;
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
title('Covariance spectrum');

figure;
plot(X(:,1:min(10,M)), 'LineWidth', 1);
grid on;
xlabel('Symbol index');
ylabel('Node value');
title('First few virtual-node traces across symbols');

%% Optional: node offsets within one symbol
tau = (0:M-1) * dt_expected;

% save('node_correlation_results_paper_consistent.mat', ...
%     'X','R','C','corrVsSep','eigVals','participationRatio', ...
%     'effectiveRank','tau','T','M','dt');

%% 
Xc = X - mean(X,1);
C = (Xc' * Xc) / size(Xc,1);

eigVals = eig(C);
eigVals = sort(real(eigVals), 'descend');

% Participation-ratio effective dimension
Deff = (sum(eigVals)^2) / sum(eigVals.^2);

% Effective rank
p = eigVals / sum(eigVals);
reff = exp(-sum(p .* log(p + eps)));

fprintf('Effective dimension (participation ratio): %.4f\n', Deff);
fprintf('Effective rank: %.4f\n', reff);