function [t_common, B_all] = runOneCase( ...
    caseInfo, ...
    N_i, input_time_points, ...
    smoldyn_visualization_on, smoldyn_flag, smoldynCmd, ...
    smold_only_data, datadirname, ...
    variableNo, variableValues, ...
    KON, KOFF, DIFFMESSENGER, DIFFRECEPTOR, DIFFRECEPTORACT, ...
    START_TIME, STOP_TIME, TIME_STEP, BIT_INTERVAL, numRECEPTOR, ...
    TXRXDISTANCE, BOUNDARYLENGTH, TXPOSITION, RXRADIUS, ...
    RXRECEPTIONSPACETHICKNESS, SAMPLING_PERIOD, N, numIterations)

    fprintf('\nRunning case: %s, r_tx = %.2f\n', caseInfo.txMode, caseInfo.r_tx)

    if smold_only_data == 1
        dirname = sprintf('%s/', datadirname);
    else
        dirname = caseInfo.dirname;
        mkdir(dirname)

        generateSimConfig_MCConfirmation(N_i, input_time_points, ...
            caseInfo.txMode, caseInfo.r_tx, smoldyn_visualization_on)

        status = copyfile('MCpointTXsphRX_generated.txt', dirname);
        if ~status
            error('Smoldyn execution file could not be copied to new directory.')
        end
    end

    %% Run Smoldyn simulations, same structure as old working code

    for j = 1:numel(variableValues)

        for i = 1:numIterations

            command = sprintf([smoldynCmd ' %sMCpointTXsphRX_generated.txt ' smoldyn_flag ' ', ...
                '--define index=%d ', ...
                '--define variableNo=%d ', ...
                '--define variableValue=%.12g ', ...
                '--define KON=%.12g ', ...
                '--define KOFF=%.12g ', ...
                '--define DIFFMESSENGER=%.12g ', ...
                '--define DIFFRECEPTOR=%.12g ', ...
                '--define DIFFRECEPTORACT=%.12g ', ...
                '--define START_TIME=%.12g ', ...
                '--define STOP_TIME=%.12g ', ...
                '--define TIME_STEP=%.12g ', ...
                '--define BIT_INTERVAL=%.12g ', ...
                '--define numRECEPTOR=%d ', ...
                '--define TXRXDISTANCE=%.12g ', ...
                '--define BOUNDARYLENGTH=%.12g ', ...
                '--define TXPOSITION=%.12g ', ...
                '--define RXRADIUS=%.12g ', ...
                '--define RXRECEPTIONSPACETHICKNESS=%.12g ', ...
                '--define SAMPLING_PERIOD=%d '], ...
                dirname, i, variableNo, variableValues(j), KON, KOFF, DIFFMESSENGER, ...
                DIFFRECEPTOR, DIFFRECEPTORACT, START_TIME, STOP_TIME, TIME_STEP, BIT_INTERVAL, numRECEPTOR, ...
                TXRXDISTANCE, BOUNDARYLENGTH, TXPOSITION, RXRADIUS, RXRECEPTIONSPACETHICKNESS, SAMPLING_PERIOD);

            fprintf('Trial %d/%d\n', i, numIterations)
            fprintf('%s\n', command)

            if smold_only_data ~= 1
                [status, cmdout] = system(command);
                if status ~= 0
                    fprintf('Smoldyn returned nonzero status = %d\n', status)
                    fprintf('%s\n', cmdout)
                end
            end
        end
    end

    %% Read simulation results, same style as old code

    ReceptorActives = cell(numel(variableValues), numIterations);
    times2 = cell(numel(variableValues), numIterations);

    t_common = [];
    B_all = [];

    for j = 1:numel(variableValues)

        for k = 1:numIterations

            fprintf('Reading varValue = %g -- numIter = %d\n', variableValues(j), k)

            allMoleculesFilename = sprintf('%sallmolecules_varNo_%d_varValue_%g_iter_%d.txt', ...
                dirname, variableNo, variableValues(j), k);

            if ~isfile(allMoleculesFilename)
                allMoleculesFilename_alt = sprintf('%sallmolecules_varNo_%d_varValue_%d_iter_%d.txt', ...
                    dirname, variableNo, round(variableValues(j)), k);

                if isfile(allMoleculesFilename_alt)
                    allMoleculesFilename = allMoleculesFilename_alt;
                else
                    error('Missing Smoldyn output file: %s', allMoleculesFilename)
                end
            end

            data2_temp = importdata(allMoleculesFilename, ' ', 1);
            data2 = data2_temp.data;

            times2{j,k} = data2(:,1);
            ReceptorActives{j,k} = data2(:,end);

            t_trial = times2{j,k};
            b_trial = ReceptorActives{j,k} / N;

            if isempty(t_common)
                t_common = t_trial;
                B_all = nan(length(t_common), numIterations);
            end

            if length(t_trial) ~= length(t_common) || max(abs(t_trial - t_common)) > 1e-12
                b_trial = interp1(t_trial, b_trial, t_common, 'linear', 'extrap');
            end

            B_all(:,k) = b_trial;
        end
    end
end