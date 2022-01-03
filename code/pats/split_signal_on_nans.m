function parts = split_signal_on_nans(sig_matrix)
    [m, n] = size(sig_matrix);
    if m < n
        sig_matrix = sig_matrix';
    end
    
    % Extract non-nan phases from signal
    nan_state_change = diff(any(isnan(sig_matrix), 2));
    nan_phase_end = find(nan_state_change == -1);
    nan_phase_start = find(nan_state_change == 1);
    nan_phase_start = [nan_phase_start; length(sig_matrix)-1];
    if ~any(isnan(sig_matrix(1, :)), 2)
        nan_phase_end = [0; nan_phase_end];
    end
    parts = cell(1, length(nan_phase_end));
    for ii = 1:length(nan_phase_end)
        parts{ii} = sig_matrix(nan_phase_end(ii)+1:nan_phase_start(ii), :);
        assert(sum(isnan(parts{ii}(:))) == 0)
    end
    % omit phases that are too short to be analyzed
    parts = parts(cellfun(@length, parts) >= 10);
end
