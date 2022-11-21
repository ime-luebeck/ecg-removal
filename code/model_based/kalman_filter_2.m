function [removedKalmanFilter, removedKalmanSmoother, phase] = ...
    kalman_filter_2(signal, rWave, fs, varargin)
%KALMANFILTER Model-based ECG-removal algorithm. This function estimates a
%noise free ECG-signal out of a cardiac contaminated EMG-signal and
%subtracts it from the measurement to clean the EMG-signal.
%   INPUT:  
%           signal  ->  cardiac contaminated EMG signal
%           rWave   ->  detected R-peaks, vector of the length of "signal"
%                       with 1 -> detected R-peak and 0-> else (see function "PeakDetection2")
%           fs      ->  Sampling rate of both signal and rWave
%           m       ->  maximum number of Gaussian kernels that are used to
%                       model the ECG-waveform (optional, default 13)
%           tau     ->  Kalman filter forgetting time. tau=[] for no
%                       forgetting factor (default)
%           gamma   ->  observation covariance adaptation-rate. 0<gamma<1
%                       and gamma=1 for no adaptation. Optional, default: 0.6.
%   OUTPUT:
%           remvodeKalmanFilter    ->  denoised ECG is estimated by the Extended
%                                      Kalman Filter and subtracted from "signal"
%           removedKalmanSmoother  ->  denoised ECG is estimated by the Extended 
%                                      Kalman Smoother and subtracted from "signal"
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

    if nargin > 3
        m = varargin{1};
    else
        m = 13;
    end

	if nargin > 4
		tau = varargin{2};
	else
		tau = [];                   % Kalman filter forgetting time. tau=[] for no forgetting factor
	end
	
	if nargin > 5
		gamma = varargin{3};
    else
        % I am unsure right now whether this should be tuned depending on fs?
		gamma = 0.6;                  % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
	end
    
    bins = round(fs/3);
    
    rWave = rWave(:)';
    [phase, ~] = PhaseCalculation(rWave);
    rWaveIndex = find(rWave);
    
    [ECGmean, ECGsd, meanphase] = MeanECGExtraction(signal, phase, bins, 1);
    [OptimumParams, ~, N] = ecg_beat_fitter_2(ECGmean, meanphase, m);
    
    fm = fs ./ diff(rWaveIndex);          % heart-rate
    w = mean(2*pi*fm);                  % average heart-rate in rads.
    wsd = std(2*pi*fm, 1);               % heart-rate standard deviation in rads.
    
    signal = signal(:)';
    
	y = [phase ; signal];
	
	X0 = [-pi 0]';              
	P0 = [(2*pi)^2 0 ;0 (10*max(abs(signal))).^2];   % Init state covariance
    
	Q = diag([(.04*OptimumParams(1:N)).^2, ...
			  (.02*ones(1,N)).^2, ...
			  (.02*ones(1,N)).^2, ...
			  (wsd)^2, ...
			  (.1*mean(ECGsd(1:round(bins/10))))^2
			  ]);
	R = [(w/fs).^2/12, 0 ; ...
		0, 0.3*(mean(ECGsd(1:round(bins/10)))).^2];

    Wmean = [OptimumParams, w, 0]';
    Vmean = [0, 0]';
    Inits = [OptimumParams, w, fs];

    InovWlen = ceil(0.2 * fs);     % innovations monitoring window length
	
    RadaptWlen = ceil(0.1 * fs);    % window length for observation covariance adaptation

	[Xekf, Xeks] = kalman_smoother_2(y, X0, P0, Q, R, Wmean, Vmean, Inits, InovWlen, tau, gamma, RadaptWlen, 1);
    
    phase = Xekf(1,:);
    Xekf = Xekf(2,:);
    Xeks = Xeks(2,:);

    removedKalmanFilter = signal - Xekf;
    removedKalmanFilter = removedKalmanFilter(:);
    
    removedKalmanSmoother = signal - Xeks;
    removedKalmanSmoother = removedKalmanSmoother(:);
end

