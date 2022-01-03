function plot_frequency_spectrum(fs, filename, varargin)
%PLOT_FREQUENCY_SPECTRUM 
%   Plots frequency spectra of multiple input signals, calculated by
%   Welch's method.
%   INPUT:  fs      ->  sampling frequency e.g. 1024 Hz
%           filename->  name of file to export plot to
%           data    ->  nx5 cell-matrix, where n is the no. of signals. The
%                       first two columns contain 2 channels of measured
%                       EMG signals and detected r peaks. Column 3 contains the
%                       measured airway pressur, 4 the time and 5 the
%                       subject number.
%           name    ->  String of name for plot legend.
%   example:    calcFrequencySpectrumOnMeasuredSignals(1024, ...
%                   dataWithRPeaks, 'Measured', ... 
%                   dataTemplatesubtraction, 'TS')
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

if mod(nargin, 2) ~= 0 || nargin < 4
   error('Invalid number of input arguments')
end
n = size(varargin{1}, 1);
% Loop over subjects
for j = 1:n
    LegendInfo = cell(1, (nargin - 2) / 2);
    fig = figure('name', ['Frequency spectrum ', 'Subject: ', varargin{1}{j,5}]);
    ax(1) = subplot(2, 1, 1);
    title(ax(1), 'Channel 1')
    ax(2) = subplot(2, 1, 2);
    title(ax(2), 'Channel 2')
    % Loop over algorithms
    for i = 1:2:nargin - 2
        LegendInfo{ceil(i / 2)} = varargin{i + 1};
        signal = varargin{i};
        channel1 = signal{j, 1}(:, 1);
        channel2 = signal{j, 2}(:, 1);
        subplot(ax(1))
        hold on
        [spectrum, f] = calcFrequencySpectrum(channel1, fs, 100, 1024);
        plot(f(2:end), 20*log10(spectrum(2:end)), 'LineWidth', 1);
        axis tight
        box on
        xlabel('Frequency (Hz)')
        ylabel('PSD (dB/Hz)')
        subplot(ax(2))
        hold on
        [spectrum, f] = calcFrequencySpectrum(channel2, fs, 100, 1024);
        plot(f(2:end), 20*log10(spectrum(2:end)), 'LineWidth', 1);
        axis tight
        box on
        xlabel('Frequency (Hz)')
        ylabel('PSD (dB/Hz)')        
    end
    subplot(ax(1))
    legend(LegendInfo)
    print_fig_to_png(fig, [filename, '_subject', varargin{1}{j,5}], 10, 10);
end

end

