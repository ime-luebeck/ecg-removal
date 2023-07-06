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

close all
clear

addpath('analysis')
addpath('ecg_utils')
addpath('emd')
addpath('filters')
addpath('format_specific_utils')
addpath('model_based')
addpath('pats')
addpath('swt_denoising')
addpath('template_subtraction')

fprintf('Load data... \n')
if isempty(dir('../data/csv/subject1.csv')) || isempty(dir('../data/csv/subject2.csv'))
    if ~exist('../data/csv', 'dir')
        mkdir('../data/csv');
    end
    unzip('../data/subject1.zip','../data/csv')
    unzip('../data/subject2.zip','../data/csv')
end
% The following call will load a total of four EMG signals, two each from
% two subjects, into a cell array of the following structure:
%       data =      {channel1, channel2, pressure, time,'01'     -> of testperson 01
%                    channel1, channel2, ...     ,     ,'02'}    -> of testperson ...
% Note that this specific structure is not necessary to run the algorithms
% in this toolbox; it is only necessary if one wants to use the
% 'filter_signals' utility function for handling multiple signals at the
% same time. You can also simply apply the removal functions to individual
% signals; refer to the readme for an example of high to do this.
rawData = load_measured_signals('../data/csv');

fs = 1024;  % signal sampling rate
fpl = 50;  % powerline interference frequency; will be taken into account / suppressed during various filtering steps
use_filtfilt = true;

%% Preprocessing and R-Peak Detection
fprintf('\n Preprocessing... \n')
% Power line removal
algoPL = @(sig) butter_filt_stabilized(sig, [fpl-1 fpl+1], fs, 'stop', use_filtfilt, 2);
tempData = filter_signals(rawData, algoPL, [1 2]);
% Low-pass against aliasing
algoLP = @(sig) butter_filt_stabilized(sig, 511.9999, fs, 'low', use_filtfilt, 6);
tempData = filter_signals(tempData, algoLP, [1 2]);
% Smooth pressure signal
algoAF = @(sig) movmean(sig, 5, 'omitnan');
tempData = filter_signals(tempData, algoAF, 3);
fprintf('\n R-Peak Detection... \n')
algoRPeak = @(sig, time) [sig, peak_detection(butter_filt_stabilized(sig, 5, fs, 'high', use_filtfilt, 6), fs, time)'];
dataWithRPeaks = filter_signals(tempData, algoRPeak, [1 2], 0, 4);

fprintf('\n Baseline Removal... \n')
% 2nd order 10 Hz high-pass for BL removal + some ECG removal
hp10 =  @(sig, rpeaks) [butter_filt_stabilized(sig, 10, fs, 'high', use_filtfilt, 2), rpeaks];
dataWithRPeaksHP10 = filter_signals(dataWithRPeaks, hp10, [1 2], 2);

% 6th order 20 Hz high-pass for BL removal + some ECG removal
hp20 =  @(sig, rpeaks) [butter_filt_stabilized(sig, 20, fs, 'high', use_filtfilt, 6), rpeaks];
dataWithRPeaksHP20 = filter_signals(dataWithRPeaks, hp20, [1 2], 2);

%% ECG-Artifact Removal
% Template Subtraction (TS)
fprintf('\n Template Subtraction (TS)... \n')
name = fullfile('../results', 'dataTS.mat');
if isempty(dir(name))
    dataTS = filter_signals(dataWithRPeaksHP20, ...
        @template_subtraction, [1 2], 2);
    save(name, 'dataTS');
else
    load(name)
end

% Adaptive Template Subtraction (ATS)
fprintf('\n Adaptive Template Subtraction (ATS)... \n')
name = fullfile('../results', 'dataATS.mat');
if isempty(dir(name))
    ats = @(signal, rpeaks) adaptive_template_subtraction(signal, rpeaks, fs);
    dataATS = filter_signals(dataWithRPeaksHP20, ats, [1 2], 2);
    save(name, 'dataATS');
else
    load(name)
end

% High-pass filter with 200 Hz cutoff frequency (HP200)
fprintf('\n 200 Hz high-pass filter (HP200)... \n')
name = fullfile('../results', 'dataHP200.mat');
if isempty(dir(name))
    dataHP200 = ...
		filter_signals(dataWithRPeaks, ...
            @(signal) butter_filt_stabilized(signal, 200, fs, 'high', use_filtfilt, 3), [1 2]);
    save(name, 'dataHP200');
else
    load(name)
end

% Adaptive Template Subtraction and subsequent high-pass filtering (ATS + HP200)
name = fullfile('../results', 'dataATS_HP200Last.mat');
if isempty(dir(name))
    dataATS_HP200Last = filter_signals(dataATS, ...
        @(signal) butter_filt_stabilized(signal, 200, fs, 'high', use_filtfilt, 3), [1 2]);
    save(name, 'dataATS_HP200Last');
else
    load(name)
end

% Stationary Wavelet Transform (SWT)
fprintf('\n Wavelet Transform (SWT)... \n')
name = fullfile('../results', 'dataSWT.mat');
if isempty(dir(name))
    swt = @(signal, rPeaks) ...
		swtden(signal, rPeaks, fs, 'h', 3, 'db2', 4.5);
    dataSWT = filter_signals(dataWithRPeaks, swt, [1 2], 2);
    save(name, 'dataSWT');
else
    load(name)
end

% Second order Extended Kalman Smoother (EKS2)
fprintf('\n Extended Kalman Smoother (EKS2)... \n')
name = fullfile('../results', 'dataEKS2.mat');
if isempty(dir(name))
    ekf2 = @(signal, rpeaks) kalman_filter_2(signal, rpeaks, fs);
    [~, dataEKS2] = ...
		filter_signals(dataWithRPeaksHP10, ekf2, [1 2], 2);
    save(name, 'dataEKS2');
else
    load(name)
end

% 25th order Extended Kalman Smoother (EKS25)
fprintf('\n Extended Kalman Smoother (EKS25)... \n')
name = fullfile('../results', 'dataEKS25.mat');
if isempty(dir(name))
    ekf25 = @(signal, rpeaks) kalman_filter_25(signal, rpeaks, fs);
    [~, dataEKS25] = ...
		filter_signals(dataWithRPeaksHP10, ekf25, [1 2], 2);
    save(name, 'dataEKS25');
else
    load(name)
end

% Empirical Mode Decomposition (EMD)
fprintf('\n Empirical Mode Decomposition (EMD)... \n')
name = fullfile('../results', 'dataEMD.mat');
if isempty(dir(name))
    emd = @(signal) EMD(signal, 60, 1);
    dataEMD = filter_signals(dataWithRPeaksHP20, emd, [1 2]);
    save(name, 'dataEMD');
else
    load(name)
end

% Probabilistic Adaptive Template Subtraction (PATS)
fprintf('\n Probabilistic Adaptive Template Subtraction (PATS)... \n')
name = fullfile('../results', 'dataPATS.mat');
if isempty(dir(name))
    pats = @(signal, rpeaks) PATS(signal, rpeaks, fs, 0, [], [], fpl);
    [~, dataPATS] = ...
        filter_signals(dataWithRPeaksHP20, pats, [1 2], 2);
    save(name, 'dataPATS');
else
    load(name)
end


%% Results
plot_results(dataWithRPeaks, [], '../results/raw_data')

fprintf('\n Plot results... \n')
plot_results(dataWithRPeaks, [], '../results/denoised_data', ...
    'TS', dataTS, ...
    'ATS', dataATS, ...
    'HP200', dataHP200, ...
    'ATS+HP200', dataATS_HP200Last, ...
    'EKS2', dataEKS2, ...
    'EKS25', dataEKS25, ...
    'SWT', dataSWT,...
    'EMD', dataEMD, ...
	'PATS', dataPATS);


%% Analysis
fprintf('\n Evaluate results...\n')
plot_pm_snr_improvement(dataWithRPeaksHP20, fs, '../results/pm_snr', ...
    dataTS, 'TS', ...
    dataATS, 'ATS', ...
    dataHP200, 'HP200', ...
    dataATS_HP200Last, 'ATS+HP200', ...
    dataEKS2, 'EKS2', ...
    dataEKS25, 'EKS25', ...
    dataSWT, 'SWT', ...
    dataEMD, 'EMD', ...
	dataPATS, 'PATS');

plot_frequency_spectrum(fs, '../results/spectrum', ...
    dataWithRPeaks, 'raw', ...
    dataTS, 'TS', ...
    dataATS, 'ATS', ...
    dataHP200, 'HP200', ...
    dataATS_HP200Last, 'ATS+HP200', ...
    dataEKS2, 'EKS2', ...
    dataEKS25, 'EKS25', ...
    dataSWT, 'SWT', ...
    dataEMD, 'EMD', ...
	dataPATS, 'PATS');