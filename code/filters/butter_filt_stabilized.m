function [signal_filt, sig_win] = butter_filt_stabilized(signal, freq_hz, fs, ftype, use_filtfilt, order)
%BUTTER_FILT_STABILIZED Stabilized Butterworth filter implementation
%   Second order butterworth filter.
%   freq_hz  ->  cut off frequency or frequency band of interest in Hz
%   fs  ->  sampling frequency
%   ftype -> filter type, e.g., "stop", "pass", "low", "high".
%   use_filtfilt -> do forward-backward filtering (doubling effective filter order!) or not (default:true)
%   order -> filter order, default is 9
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

    assert(~any(isnan(signal)), 'Does not work with nan since it is an IIR filter. Consider using fillmissing or an FIR filter.');

    if size(signal, 1) > size(signal, 2)
        signal = signal';
        was_col_vec = true;
    else
        was_col_vec = false;
    end
    
    if nargin < 5
        use_filtfilt = true;
    end
    
    if nargin < 6
        order = 9;
    end
    
    % Filter produces initialization artifacts at beginning & end of data set.
    % To reduce these, apply a windowing function that creates a smooth transition to/from zero.
    % Choose the length of the initialization section adaptively, such that it covers two oscillations
    % of the lower cutoff frequency (up to a maximum of 120s or half of the signal; whichever is smaller).
    half_win_len = min(round(2 * fs / min(freq_hz)), min(round(fs * 120), floor(length(signal)/2)));
    hann_win = hann(half_win_len*2)';
    sig_win = [hann_win(1:half_win_len), ones(1, length(signal) - 2*half_win_len), hann_win(half_win_len+1:end)];
    signal_windowed = sig_win .* signal;
    
    % Requires Signal Processing Toolbox
    [A, B, C, D] = butter(order, freq_hz / (fs/2), ftype);
    
    % Plot bode diagram of that filter (probably there's also a simpler way?)
    % [b, a] = ss2tf(A, B, C, D);
    % bode(tf(b, a, 1/fs))
    
    % Standard SS form is numerically instable, hence transform to more stable
    % but equivalent second-order section form.
    [sos, g] = ss2sos(A, B, C, D);
    
    % The following call to filtfilt always produces a warning because of
    % ambiguous arguments, but it does exactly what it's meant to do.
    w = warning ('off','all');    
    
    if use_filtfilt
        % Perform zero-phase filtering [doubles effective filter order!]
        signal_filt = filtfilt(sos, g, signal_windowed);
    else
        % Perform standard forward filtering
        signal_filt = filter(sos, g, signal_windowed);
    end
    
    % Restore original warning state
    warning(w);    
    
    if was_col_vec
        signal_filt = signal_filt';
    end
end