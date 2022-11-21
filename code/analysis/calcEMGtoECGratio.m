function SNR = calcEMGtoECGratio(dataWithRPeaks, fs, varargin)
%CALCEMGTOECGRATIO Calculates the Signal-to-Noise ratio of cardiac
%contaminated EMG signals
%
%   The SNR is calculated by the RMS of the active EMG parts divided by the
%   the RMS of QRS-complexes in the signal. Active parts in the EMG signal
%   are determined by a simultaneously measured pressure signal
%   (inspiration -> active).
%   INPUT:  dataWithRPeaks  ->  nx5 cell-matrix, where n is the no. of
%                               signals. The first two columns contain 2
%                               channels of measured EMG-signals and detected r peaks. Column 3
%                               contains the measured airway pressure, 4 the time
%                               and 5 the subject number.
%           fs              ->  Sampling rate of all signals
%           flagPlot        ->  If 1, results are plotted.
%           data            ->  SNR is calculated of alternative set of data with
%                               same structure as dataWithRPeaks.
%           example: calcEMGtoECGratio(dataWithRPeaks, 1024, 1, dataTemplateSubtraction);
%   OUTPUT: SNR             ->  mx2 matrix with SNR of 2 channels for m subjects 
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

if nargin > 2 && ~isempty(varargin{1})
    flagPlot = varargin{1};
else
    flagPlot = 1;
end

if nargin > 3 && ~isempty(varargin{2})
    data = varargin{2};
else
    data = dataWithRPeaks;
end

if nargin > 4 && ~isempty(varargin{3})
    skip_seconds = varargin{3};
else
    skip_seconds = 0;
end

n = size(data, 1);
N = size(data{1, 1}, 1);

skip_samples = round(fs*skip_seconds);

SNR = zeros(n, 2);
for i = 1:n
    % Skip first + last 10s as these often contain artifacts
    signal_1 = data{i, 1}(skip_samples+1:end-skip_samples, 1); % EMG channel 1
    r_peak_locs_1 = find(dataWithRPeaks{i, 1}(skip_samples+1:end-skip_samples, 2)); % Index of R peaks in channel 1
    signal_2 = data{i, 2}(skip_samples+1:end-skip_samples, 1); % EMG channel 2
    r_peak_locs_2 = find(dataWithRPeaks{i, 2}(skip_samples+1:end-skip_samples, 2)); % Index of R peaks in channel 2
    
    % Extract insp phases - we assume these to be phases of muscular activity
    pressure = movmean(dataWithRPeaks{i, 3}, round(0.1*fs), 'omitnan');
    insp_mask = detectInspFromPressure(pressure(skip_samples+1:end-skip_samples), fs);
    
    % Calculate empirical SNR
    [snr1, activity_only_bl_removed_1, activity_mask_1] = SNRmeasured(signal_1, ... 
					   r_peak_locs_1, fs, insp_mask);
    [snr2, activity_only_bl_removed_2, activity_mask_2] = SNRmeasured(signal_2, ... l
					   r_peak_locs_2, fs, insp_mask);
                   
    SNR(i, :) = [snr1, snr2];

    if flagPlot == 1
        figure;
        
        h1 = subplot(3, 2, 1);
        plot(signal_1);
        hold on;
        scale_1 = mean(abs(signal_1));
        plot(r_peak_locs_1, signal_1(r_peak_locs_1), 'ro');
        plot(scale_1*activity_mask_1, 'black');
        plot(-scale_1*insp_mask + 0.2*scale_1, 'magenta');
        
        h2 = subplot(3, 2, 2);
        plot(signal_2);
        hold on;
        scale_2 = mean(abs(signal_2));
        plot(r_peak_locs_2, signal_2(r_peak_locs_2), 'ro');
        plot(scale_2*activity_mask_2, 'black');
        plot(-scale_2*insp_mask + 0.2*scale_2, 'magenta');
        
        h3 = subplot(3, 2, 3);
        plot(find(activity_mask_1), activity_only_bl_removed_1);
        
        h4 = subplot(3, 2, 4);
        plot(find(activity_mask_2), activity_only_bl_removed_2);
        
        h5 = subplot(3, 2, 5);
        plot(pressure);
        hold on;
        pressure_scale = mean(abs(pressure - max(pressure)));
        plot(max(pressure) - pressure_scale*insp_mask, 'magenta');
        
        h6 = subplot(3, 2, 6);
        plot(pressure);
        hold on;
        plot(max(pressure) - pressure_scale*insp_mask, 'magenta');
        
        linkaxes([h1, h3, h5], 'x');
        linkaxes([h2, h4, h6], 'x');
    end
end
    
end
