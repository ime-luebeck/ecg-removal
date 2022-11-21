function [SNR, activity_only, activity_only_mask, signal_bl_removed] = ...
    SNRmeasured(signal, rWaveIndex, fs, activity_mask)
%SNRMEASURED Calculating SNR (EMG-to-ECG ratio) of a cardiac contaminated
%EMG signal (optinally) with the help of a simultaneosly measured resp activity signal.
%
%   Active EMG sections are extracted during inspiration and QRS-complexes
%   are extracted during expiration. Sections in the signal exceeding a
%   threshold are excluded from SNR calculation. This affects only
%   processed signals with outliers. The ratio is formed by RMS.
%   
%   INPUT:  signal        ->  Vector with cardiac contaminated EMG signal.
%           rWaveIndex    ->  Indices of detected r peaks in 'signal'.
%           fs            ->  Sampling rate of signal (Hz)
%           activity_mask ->  Vector with same length as 'signal' containing
%                             1 where EMG activity is present, and 0
%                             elsewhere.
%                             This is optional. If it is not provided,
%                             "Signal" is anything where there is no QRS,
%                             and "Noise" power is calculated over all QRS
%                             regions.
%   OUTPUT: SNR           ->  Scalar with EMG-to-ECG ratio of 'signal'.
%
% Copyright 2021 Institute for Electrical Engineering in Medicine, 
% University of Luebeck
% Eike Petersen, Julia Sauer
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

[~, QRS_mask] = extractQRS(signal, rWaveIndex, fs);

% Drop samples where both QRS complex and EMG are active
if nargin > 3
    to_drop_mask = QRS_mask & activity_mask;
    QRS_mask(to_drop_mask) = 0;
    activity_only_mask = activity_mask;
    activity_only_mask(to_drop_mask) = 0;
else
    % No knowledge about muscle activity provided. We can only
    % differentiate between QRS yes and no.
    activity_only_mask = ~QRS_mask;
end

QRS_only = signal(QRS_mask~=0);

% Extracted regions of EMG activity will still contain strong P wave
% influence. To suppress this, filter the complete signal with a zero-phase high-pass 
bl_filt_cutoff = 30; % Hz
use_filtfilt = 1;
order = 6;
signal_bl_removed = butter_filt_stabilized(signal, bl_filt_cutoff, fs, 'high', use_filtfilt, order);
outlier_cutoff = quantile(abs(signal_bl_removed(activity_only_mask~=0)), 0.999);
% reject 20ms window around each outlier
twenty_ms_samples = round(0.01*fs);
outlier_mask = movmax(abs(signal_bl_removed) > outlier_cutoff, twenty_ms_samples);
activity_only_mask = activity_only_mask & ~outlier_mask;
activity_only = signal_bl_removed(activity_only_mask~=0);


SNR = mean(activity_only.^2) / mean(QRS_only.^2);

end