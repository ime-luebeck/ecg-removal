function [ECGremovedSignal, ECGremovedSignal_white_plr_hp] = PATS(signal, rPeaks, fs, varargin)
% PATS Probabilistic Adaptive Template Subtraction that takes time variability and the influence of the RR
% interval on the heart beat signal into account
%
% For a detailed motivation, derivation, and description of this algorithm, refer to Chapter 5 of 
% Eike Petersen (2022), "Model-based Probabilistic Inference for Monitoring Respiratory Effort
% using the Surface Electromyogram" (Dissertation, forthcoming).
%
% Copyright 2021 Institute for Electrical Engineering in Medicine, 
% University of Luebeck
% Eike Petersen
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the 
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
% THE USE OR OTHER DEALINGS IN THE SOFTWARE.


if nargin > 3 && ~isempty(varargin{1})
    variant = varargin{1};
else
    variant = 4;
end

if nargin > 4 && ~isempty(varargin{2})
    use_poly = varargin{2};
else
    use_poly = true;
end

if nargin > 5 && ~isempty(varargin{3})
    delays = varargin{3};
else
    delays = 0;
end

% Power line frequency (Hz); will be taken into account and suppressed in some filtering steps
if nargin > 6 && ~isempty(varargin{4})
    fpl = varargin{4};
else
    fpl = 50;
end

rPeakIdces = find(rPeaks(:));

signal = signal(:);
[N, m] = size(signal);
numQRS = length(rPeakIdces);

subtractionSignal = zeros(N, 1);

% distance between R waves
distance = rPeakIdces(2:end) - rPeakIdces(1:end-1);

% center between QRS complexes (row vector)
center = ceil(distance/2) + rPeakIdces(1:end-1);
center = [1; center(:); length(signal)+1];

template_len = round(1.2 * mean(distance));
% ensure an odd template length, such that the R peak is located exactly in the middle
if mod(template_len, 2) == 0
    template_len = template_len + 1;
end
template_R_peak_loc = ceil(template_len/2);
template_half_len = floor(template_len/2);


% max template length
length_left = rPeakIdces(1:end) - center(1:end-1);
length_right = center(2:end) - 1 - rPeakIdces(1:end);
max_length_left = max(length_left);
max_length_right = max(length_right);
max_length = max_length_left + max_length_right + 1;
beat_mat_r_peak_idx = max_length_left + 1;

if variant > 1
    % Construct ECG beat matrix
    beats = zeros(numQRS-2, max_length);
    for ii = 1:numQRS
        ldiff = rPeakIdces(ii)-max_length_left - 1;
        rdiff = N - (rPeakIdces(ii)+max_length_right);
        if ldiff < 0
            beats(ii, 1-ldiff:end) = signal(1:rPeakIdces(ii)+max_length_right)';
        elseif rdiff < 0
            beats(ii, 1:end+rdiff) = signal(rPeakIdces(ii)-max_length_left:N)';
        else
            beats(ii, :) = signal(rPeakIdces(ii)-max_length_left:rPeakIdces(ii)+max_length_right)';
        end
    end
    mean_beat = mean(beats);    
end

if variant == 0
    % No external shape information input signal
    input_signals = [];
    n_features = 0;
    
else
    % Construct shape information signal
    
    if variant == 1
        % Use RR interval
        shape_sig_raw = [distance(1); distance];

    elseif variant == 2
        % Use ECG amplitude
        shape_sig_raw = max(abs(beats), [], 2); % row maxima

    elseif variant == 3
        % Use both: RR interval and ECG amplitude
        dist_sig = [distance(1); distance];
        amp_sig = max(abs(beats), [], 2);
        shape_sig_raw = [dist_sig, amp_sig];

    elseif variant == 4
        % Use first 3 principal components
        qrs_half_width = round(0.1*fs); % Take 100ms left and right
        qrs_complexes = beats(:, beat_mat_r_peak_idx-qrs_half_width:beat_mat_r_peak_idx+qrs_half_width);
        % Ignore first + last five beats as these often contain initialization artifacts
        [~, ~, V] = svd(qrs_complexes(6:end-5, :), 'econ');
        % Calculate projection of (the QRS complex of) each beat onto first three
        % PCs and use these as the shape parameters
        shape_sig_raw = qrs_complexes * V(:, 1:3);

    elseif variant == 5
        % Use first 3 principal components and RR interval
        qrs_half_width = round(0.1*fs); % Take 100ms left and right
        qrs_complexes = beats(:, beat_mat_r_peak_idx-qrs_half_width:beat_mat_r_peak_idx+qrs_half_width);
        % Ignore first + last five beats as these often contain initialization artifacts
        [~, ~, V] = svd(qrs_complexes(6:end-5, :), 'econ');
        % Calculate projection of (the QRS complex of) each beat onto first three
        % PCs and use these as the shape parameters
        shape_sig_raw = [qrs_complexes * V(:, 1:3), [distance(1); distance]];

    else
        error('Unknown variant requested.');
    end

    % center and normalize
    shape_sig = shape_sig_raw - mean(shape_sig_raw(6:end-5, :));
    shape_sig = shape_sig ./ std(shape_sig(6:end-5, :));


    if use_poly
        if size(shape_sig_raw, 2) == 1
            features = @(shape) [shape.^5, shape.^4, shape.^3, shape.^2, shape, ones(size(shape, 1), 1)];
            n_features = length(features(0));
        elseif size(shape_sig_raw, 2) == 2
            features = @(shape) [prod(shape, 2), shape(:, 1), shape(:, 2), ones(size(shape, 1), 1)];
            n_features = length(features([0, 0]));
        elseif size(shape_sig_raw, 2) == 3
            features = @(shape) [prod(shape, 2), prod(shape(:, [1 2]), 2), prod(shape(:, [1 3]), 2), prod(shape(:, [2 3]), 2), ...
                shape(:,1), shape(:,2), shape(:,3), ones(size(shape, 1), 1)];
            n_features = length(features([0, 0, 0]));
        elseif size(shape_sig_raw, 2) == 4
            features = @(shape) [prod(shape, 2), prod(shape(:, [1 2 3]), 2), prod(shape(:, [1 2 4]), 2), prod(shape(:, [1 3 4]), 2), ...
                prod(shape(:, [1 2]), 2), prod(shape(:, [1 3]), 2), prod(shape(:, [1 4]), 2), shape(:, 1), ...
                prod(shape(:, [2 3 4]), 2), prod(shape(:, [2 3]), 2), prod(shape(:, [2 4]), 2), prod(shape(:, [3 4]), 2), ...
                shape(:, 2), shape(:, 3), shape(:, 4), ones(size(shape, 1), 1)];
            n_features = length(features([0, 0, 0, 0]));    
        else
            error('Unsupported number of dimensions requested.');
        end
    else
        features = @(shape) shape;
        n_features = size(shape_sig, 2);  
    end

    % Set up input signals
    nin = size(shape_sig, 2);
    input_signals = zeros(size(shape_sig, 1), n_features + delays*nin);
    input_signals(:, 1:n_features) = features(shape_sig);
    for ii = 1:delays
        input_signals(:, (1:nin)*ii+n_features) = [repmat(shape_sig(1, :), ii, 1); shape_sig(1:end-ii, :)];
    end
    n_features = size(input_signals, 2);  % now including delays
end


%% 1. Whiten signal (to justify assuming no correlation between subsequent measurement noise samples)
% Important: whiten using the spectrum of EMG/no-ECG sections, because those represent the
% measurement noise!
% Estimation is a bit complicated: since we only use the EMG/no-ECG sections, we have a lot of gaps
% in our signal, which standard spectral estimators cannot deal with. However, basic AR model
% identification (using prediction error minimization) in the SysID toolbox does the same thing as,
% e.g., arburg, *and* can deal with multiple "experiments" for the same identification task.
% Use that!
[~, ~, activity_only_mask, signal_bl_removed] = SNRmeasured(signal, rPeakIdces, fs);
activity_only_signal = signal_bl_removed;
activity_only_signal(~activity_only_mask) = nan;
% ignore first + last 5s as these often contain artifacts
activity_only_signal(1:round(5*fs)) = nan;
activity_only_signal(end-round(5*fs):end) = nan;
ar_id_data = non_nan_phases_iddata(activity_only_signal, [], [], fs);
sys = arx(ar_id_data, 6);
a = sys.a;
%a = arburg(activity_only_bl_removed, 6);

% filter with inverse spectrum
signal_white = filter(a, 1, signal);


%% 2. Estimate time-varying measurement (=EMG) noise covariance, using segments where there is likely no cardiac activity
[~, ~, activity_only_mask, signal_bl_removed] = SNRmeasured(signal_white, rPeakIdces, fs);
% Ignore first + last 5s as these often contain artifacts
% Identify a MOVING variance as the variance of the EMG signal depends on the level of activity
not_first_last_5s_mask = (((1:length(signal_white)) >= round(fs*5)) & ((1:length(signal_white)) <= length(signal_white) - round(fs*5)))';
emg_variance = nan * ones(size(signal_white));
emg_variance(activity_only_mask~=0 & not_first_last_5s_mask) = ...
    movvar(signal_bl_removed(activity_only_mask~=0 & not_first_last_5s_mask), round(0.3*fs));
emg_variance = fillmissing(emg_variance, 'linear', 'EndValues', 'nearest');
%emg_variance = var(activity_only_bl_removed(round(fs*5):end-round(fs*5)));


%% 3. Identify heart beat template dynamical models + time courses 

% This will contain the estimated time course of the different cardiac template samples
template_est = zeros(numQRS, template_len);

% Loop over indices of the heart beat template (roughly 1...1000)
p0 = 0.1 * ones(5+n_features, 1);

for ii = 1:template_len
    signal_idces = rPeakIdces + (ii - template_R_peak_loc);
    signal_idces = signal_idces(signal_idces > 0 & signal_idces <= N);
    
    meas_signal = signal_white(signal_idces);
    emg_var_loc = emg_variance(signal_idces);
    if rPeakIdces(1) + ii - template_R_peak_loc <= 0
        meas_signal = [nan; meas_signal];
        emg_var_loc = [nan; emg_var_loc];
    end
    
    if rPeakIdces(end) + ii - template_R_peak_loc > N
        meas_signal = [meas_signal; nan];
        emg_var_loc = [emg_var_loc; nan];
    end
    
    % 3a. Perform ML estimation to identify the state-space model for this heart beat template sample
    % Reuse the solution of the previous heartbeat for initialization
    % Replace fmincon by NOMAD once available on unix?
    % Ignore first + last five beats    
    lb = [-2, -2, -2, -0.2 * ones(1, n_features), -0.5, 1e-8]';
    ub = [2, 2, 2, 0.2 * ones(1, n_features), 0.5, 0.2]';
    opt_params = fmincon(@(par) sample_time_course_est(meas_signal(6:end-5), input_signals(6:end-5, :), par, ...
                         emg_var_loc(6:end-5)), p0, [], [], [], [], lb, ub, [], ...
                         optimoptions(@fmincon, 'Display', 'iter', 'Diagnostics', 'on', 'UseParallel', true));
%     if ispc
%         % Use OPTI toolbox version of NOMAD solver
%         opts = optiset('solver','nomad', 'display', 'iter', 'maxfeval', 1000, 'tolrfun', 1e-7, 'tolafun', 1e-4);
%         Opt = opti('fun', @(par) sample_time_course_est(meas_signal, input_signals, par, emg_variance), ...
%             'bounds', lb, ub, 'Options', opts);
%         opt_params = solve(Opt, p0);
%     elseif isunix
%         % Use self-compiled NOMAD mex file
%         opts = nomadset('display_degree', 2, 'max_bb_eval', 1000);
%         opt_params = nomad(@(par) sample_time_course_est(meas_signal, input_signals, par, emg_variance), p0, lb, ub, opts);
%     else
%         error('Unsupported OS.');
%     end
    
    % Use the optimal solution for this sample as the starting point for the next sample
    p0 = opt_params;
    
    % 3b. Use the ML-optimal model identified above to estimate the time course of this sample
    % (now for all beats, also the five first+last ones)
    [~, template_est(:, ii)] = sample_time_course_est(meas_signal, input_signals, opt_params, emg_var_loc);
end


%% 4. Combine the identified time courses to an estimated cardiac signal, taking care of gaps + overlap

for ii = 1:numQRS
    template_loc = template_est(ii, :)';
    
    % 4A. Deal with the first half of the beat, including the R peak
    if ii > 1
        overlap_left = template_len - (rPeakIdces(ii) - rPeakIdces(ii-1));
    else
        overlap_left = template_half_len - rPeakIdces(1) + 1;
    end
    
    if overlap_left < 0
        % This is an exceptionally long RR interval; there is a gap to be filled.
        % Interpolate linearly between this template's start and the previous template's end.
        subtractionSignal(rPeakIdces(ii)-template_half_len:rPeakIdces(ii)) = template_loc(1:template_R_peak_loc);
        % The component belonging to the template of the previous beat has already been added.
        stepsize = 1/(abs(overlap_left)+1);
        weights = stepsize:stepsize:(1-stepsize);
        idces_gap = rPeakIdces(ii)-template_half_len+overlap_left:rPeakIdces(ii)-template_half_len-1;
        invalid_idces = idces_gap <= 0 | idces_gap > N;        
        idces_gap = idces_gap(~invalid_idces);
        weight_idces = 1:length(weights);
        weight_idces = weight_idces(~invalid_idces);           
        subtractionSignal(idces_gap) = subtractionSignal(idces_gap) + weights(weight_idces)' * template_loc(1);
        
    elseif overlap_left > 0
        % There is some overlap between this and the previous template.
        % In the overlap region, gradually transition from one template to the other.
        % if overlap_left > template_half_len, we're in serious trouble
        % anyway... (two peaks are extremely close together).        
        overlap_left = min(overlap_left, template_half_len - 50);
        weights = ones(1, template_half_len+1);
        stepsize = 1/(overlap_left+1);
        weights(1:overlap_left) = stepsize:stepsize:(1-stepsize);
        idces_overlap = rPeakIdces(ii)-template_half_len:rPeakIdces(ii);
        invalid_idces = idces_overlap <= 0 | idces_overlap > N;
        idces_overlap = idces_overlap(~invalid_idces);
        template_idces = 1:template_R_peak_loc;
        template_idces = template_idces(~invalid_idces);
        subtractionSignal(idces_overlap) = subtractionSignal(idces_overlap) + weights(template_idces)' .* template_loc(template_idces);
        
    else
        % The very rare case where the RR interval matches the template length exactly.
        subtractionSignal(center(ii):rPeakIdces(ii)) = template_loc(1:template_R_peak_loc);
    end
    
    % 4B. Deal with the second half of the beat.
    if ii < numQRS
        overlap_right = template_len - (rPeakIdces(ii+1) - rPeakIdces(ii));
    else
        overlap_right = template_half_len - (N - rPeakIdces(end));
    end
    
    if overlap_right < 0
        % This is an exceptionally long RR interval; there is a gap to be filled.
        % Interpolate linearly between this template's end and the next template's start.
        subtractionSignal(rPeakIdces(ii)+1:rPeakIdces(ii)+template_half_len) = template_loc(template_R_peak_loc+1:end);
        % The component belonging to the template of the next beat will be added later
        stepsize = 1/(abs(overlap_right)+1);
        weights = (1-stepsize):-stepsize:stepsize;
        idces_gap = rPeakIdces(ii)+template_half_len+1:rPeakIdces(ii)+template_half_len-overlap_right;
        invalid_idces = idces_gap <= 0 | idces_gap > N;        
        idces_gap = idces_gap(~invalid_idces);
        weight_idces = 1:length(weights);
        weight_idces = weight_idces(~invalid_idces);            
        subtractionSignal(idces_gap) = weights(weight_idces)' * template_loc(end);
        
    elseif overlap_right > 0
        % There is some overlap between this and the next template.
        % In the overlap region, gradually transition from one template to the other.
        % if overlap_right > template_half_len, we're in serious trouble anyway: two peaks are very close together.
        overlap_right = min(overlap_right, template_half_len - round(0.05*fs));
        weights = ones(1, template_half_len);
        stepsize = 1/(overlap_right+1);
        weights(end-overlap_right+1:end) = (1-stepsize):-stepsize:stepsize;
        idces_overlap = rPeakIdces(ii)+1:rPeakIdces(ii)+template_half_len;
        invalid_idces = idces_overlap <= 0 | idces_overlap > N;        
        idces_overlap = idces_overlap(~invalid_idces);
        template_idces = (template_R_peak_loc+1):template_len;
        template_idces = template_idces(~invalid_idces);    
        weight_idces = template_idces - template_R_peak_loc;
        subtractionSignal(idces_overlap) = weights(weight_idces)' .* template_loc(template_idces);
        
    else
        % The very rare case where the beat matches the template length exactly.
        subtractionSignal(rPeakIdces(ii)+1:center(ii+1)-1) = template_loc(template_R_peak_loc+1:end);
    end
end


%% 5. Estimate EMG as difference and un-whiten signal
% remove estimated cardiac component
ECGremovedSignal_white = signal_white - subtractionSignal;

% use the inverse of the whitening filter from step 1 to restore EMG frequency characteristics
%ECGremovedSignal = filter(1, a, ECGremovedSignal_white);
% but make sure to remove all low-frequency + 50Hz stuff before doing that, in order not to amplify noise!
use_filtfilt = true; half_order = 2;
ECGremovedSignal_white_plr = butter_filt_stabilized(ECGremovedSignal_white, [fpl-1 fpl+1], fs, 'stop', use_filtfilt, half_order);
ECGremovedSignal_white_plr_hp = butter_filt_stabilized(ECGremovedSignal_white_plr, 20, fs, 'high', use_filtfilt, half_order);
ECGremovedSignal = filter(1, a, ECGremovedSignal_white_plr_hp);

%% 6. Done!

end


function [energy, time_course_est, time_course_cov] = ...
    sample_time_course_est(meas_signal, input_signals, model_params, emg_variance)
    
    n_features = size(input_signals, 2);
    A = [model_params(1:3)'; 1 0 0; 0 1 0];
    Q = zeros(3);
    Qmean = model_params(4+n_features);    
    Q(1,1) = model_params(5+n_features)^2;

    C = [1 0 0];
    R = emg_variance;
    
    if n_features == 0
        [~, X_smooth, ~, P_smooth, ~, ~, ~, energy] = ...
            KF_RTS(meas_signal', A, C, Q, R, [], [], 'Qmean', Qmean);   
    else
        B = [model_params(4:3+n_features)'; zeros(2, n_features)];
        [~, X_smooth, ~, P_smooth, ~, ~, ~, energy] = ...
            KF_RTS(meas_signal', A, C, Q, R, [], [], 'B', B, 'u', input_signals, 'Qmean', Qmean);        
    end
    
    if nargout > 1
        time_course_est = X_smooth(1, :)';
        time_course_cov = squeeze(P_smooth(1, 1, :));
    end
end