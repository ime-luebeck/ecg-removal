function merged_nan_parts = non_nan_phases_iddata(output_signal, input_signals, output_filter, fs)

    

    if nargin < 2 || isempty(input_signals)
        sig_mat = [output_signal];
    else
        sig_mat = [output_signal, input_signals];
    end
    
    if nargin < 3
        output_filter = [];
    end
    
    if nargin < 4
        fs = 1;
    end
    
    Ts = 1/fs;
    
    % Split into non-nan parts of the signal
    signal_parts = split_signal_on_nans(sig_mat);
    % Use the different non-nan parts as different "experiments" belonging
    % to a single iddata dataset
    signal_parts_iddata = cell(1, length(signal_parts));
    for ii = 1:length(signal_parts)
        % Construct iddata object
        curr_signal_mat = signal_parts{ii};
        curr_output = curr_signal_mat(:, 1);
        if ~isempty(output_filter)
            curr_output = output_filter(curr_output);
        end
        if nargin >= 2
            curr_input = curr_signal_mat(:, 2:end);
        else
            curr_input = [];
        end
        signal_parts_iddata{ii} = iddata(curr_output, curr_input, Ts);
    end
    merged_nan_parts = merge(signal_parts_iddata{:});
end