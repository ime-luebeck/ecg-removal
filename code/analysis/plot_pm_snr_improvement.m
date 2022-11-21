function plot_pm_snr_improvement(dataWithRPeaks, fs, filename, varargin)
%PLOT_PM_SNR_IMPROVEMENT Calculates and plots the improvement of
%periodicity measurement and signal-to-noise ratio of the ECG-removed
%signals compared to raw EMG signals.
%
% Copyright 2019 Institute for Electrical Engineering in Medicine, 
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

algo_data = varargin(1:2:end);
algo_names = varargin(2:2:end);

[SNRimprovement, ~] = analyzeSNRimprovement(dataWithRPeaks, fs, algo_data{:});
[pmImprovement, ~] = analyzePeriodicityMeasurement(dataWithRPeaks, [], algo_data{:});

fig = figure;
sp_snr = subplot(2,1,1);
title(sp_snr, 'SNR improvement')
ylabel(sp_snr, '\Delta SNR [dB]')
box on
grid on
hold on
sp_pm = subplot(2,1,2);
title(sp_pm, 'Periodicity Measure improvement')
ylabel(sp_pm, '\Delta PM [dB]')
box on
grid on
hold on

boxplot(sp_snr, SNRimprovement, 'Labels', algo_names)
boxplot(sp_pm, pmImprovement, 'Labels', algo_names)

print_fig_to_png(fig, filename, 14, 10);
end

