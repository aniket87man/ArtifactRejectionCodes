% calculates max cross-correlation of multi-trial data with mean, and corresponding lags

function [optim_shifts, max_correlations] = getCorrelations(signals, avg_signal)

    N_signal = size(signals, 1);
    T_signal = size(signals, 2);
    optim_shifts = zeros(N_signal, 1);
    max_correlations = zeros(N_signal, 1);
    if T_signal ~= size(avg_signal, 2)
        if T_signal < size(avg_signal, 2)
            signals = [signals, zeros(N_signal, size(avg_signal, 2) - T_signal)];
        else
            avg_signal = [avg_signal, zeros(1, T_signal - size(avg_signal, 2))];
        end
    end

    for i = 1:N_signal
        s_i = signals(i, :);
        [corrs, lags] = xcorr(s_i, avg_signal, 'coeff');

        [max_corr, idx] = max(corrs);
        optim_lag = lags(idx);

        optim_shifts(i) = optim_lag;
        max_correlations(i) = max_corr;
        % % % % testing
        % % if c_num >= 1 && c_num <= 14
        % %     fprintf('Signal %d: Optimal lag = %d samples, Max Correlation = %.4f\n', i, optim_lag, max_corr)
        % % end
    end
    % % % % fprintf('\nCorrelation calculation successful\n\n')

end