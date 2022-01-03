function [pmImprovement, pmMeas] = analyzePeriodicityMeasurement(dataWithRPeaks, ...
    env_method, data, varargin)
%ANALYZEPERIODICITYMEASURE Calculates Periodicity Measure of multiple input data sets. 
%The values are calculated by the function periodicityMeasure().
%   INPUT:  dataWithRPeaks  ->  nx5 cell-matrix, where n is the no. of
%                               signals. The first two columns contain 2
%                               channels of measured EMG-signals and detected r peaks. Column 3
%                               contains the measured airway pressure, 4 the time
%                               and 5 the subject number.
%           data            ->  Cell-matrix with same structure as
%                               'dataWithRPeaks' whose Periodicty Measure
%                               shall be calculated.
%           env_method         ->  String to specify if PM of envelope
%                               and by which method shall be calculated and plotted
%                               If [] or omitted, PM of signals is
%                               calculated and plotted
%                               If 'rms', envelope is calculated by RMS
%                               If 'med', Median
%                               If 'abs', Absolute
%   example:    analyzePeriodicityMeasurement('PM_data.csv', dataWithRPeaks, ...
%                   dataKalmanFilter1, ...
%                   dataKalmanFilter2)
%               or
%               analyzePeriodicityMeasurement('PM_data.csv', dataWithRPeaks, 'med' ...
%                   dataKalmanFilter1, ...
%                   dataKalmanFilter2 ...
%                   )
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

    n = size(varargin, 2) + 1;
    fs = 1024;
    skip_seconds = 10;
    
    % PM of measured signals (dataWithRPeaks)
    pmMeas = calcPeriodicityMeasures(dataWithRPeaks, dataWithRPeaks, env_method, fs, skip_seconds);
    pmMeas = pmMeas(:);
    % make db values
    pmMeas = 20*log10(pmMeas);
	
    % PM of ECG-removed signals
    disp('First algo...')
    pmData = calcPeriodicityMeasures(dataWithRPeaks, data, env_method, fs, skip_seconds);
    pmData = pmData(:);
    % make db values
    pmData = 20*log10(pmData);

    % Store results for export
    pmImprovement = zeros(length(pmMeas), n);
    pmImprovement(:, 1) = pmMeas - pmData;

	for ii = 2:n
        disp('Next algo...')
		pmData = calcPeriodicityMeasures(dataWithRPeaks, varargin{ii-1}, env_method, fs, skip_seconds);
        pmData = pmData(:);
        % make db values
        pmData = 20*log10(pmData);
        pmImprovement(:, ii) = pmMeas - pmData;
	end
end


function pm = calcPeriodicityMeasures(dataWithRPeaks, data, env_method, fs, skip_seconds)
	n = size(data, 1);
	pm = zeros(n, 2);
	% Loop over recording sessions
    skip_samples = round(fs*skip_seconds);
	for i = 1:n
        % first channel
		rPeaks1 = find(dataWithRPeaks{i,1}(skip_samples+1:end-skip_samples,2));
        rWave1 = dataWithRPeaks{i,1}(skip_samples+1:end-skip_samples,2);
		rWave1 = rWave1(:)';
		[phase1, ~] = PhaseCalculation(rWave1);
        if ~isempty(env_method)
            sig = emgEnv(data{i,1}(skip_samples+1:end-skip_samples,1), 100, env_method);
        else
            sig = data{i,1}(skip_samples+1:end-skip_samples,1);
        end
		pm1 = periodicityMeasure(sig, rPeaks1, phase1);
		
        % second channel
		rPeaks2 = find(dataWithRPeaks{i,2}(skip_samples+1:end-skip_samples,2));
        rWave2 = dataWithRPeaks{i,2}(skip_samples+1:end-skip_samples,2);
		rWave2 = rWave2(:)';
		[phase2, ~] = PhaseCalculation(rWave2);
        if ~isempty(env_method)
            sig = emgEnv(data{i,2}(skip_samples+1:end-skip_samples,1), 100, env_method);
        else
            sig = data{i,2}(skip_samples+1:end-skip_samples,1);
        end
		pm2 = periodicityMeasure(sig, rPeaks2, phase2);
		pm(i,:) = [pm1, pm2];
	end
end