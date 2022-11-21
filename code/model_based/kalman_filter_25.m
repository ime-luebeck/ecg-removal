function [EKF25, EKS25] = kalman_filter_25(signal, rWave, fs)
%EKF25 ECG-removal of EMG signals with the help of a non-linear ECG model
%and the extended kalman filter and smoother
%   INPUT:  signal  ->  Vector with EMG signal
%           rWave   ->  Vector with same length as 'signal' which is 1 when
%                       an r peak is detected and 0 else
%           fs      ->  Sampling rate of signal and rWave
%   OUTPUT: EKF25   ->  Vector with same length as 'signal' with ECG-removed
%                       EMG signal by extended kalman filter
%           EKS25   ->  Vector with same length as 'signal' with
%                       ECG-removed EMG signal by extended kalman smoother
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


% linear phase calculation
phase = PhaseCalculation(rWave);
rWaveIndex = find(rWave);

% determination of mean ECG with linear phase
bins = round(fs/3);
[meanECG, ~, meanPhase] = MeanECGExtraction(signal, phase, bins, 1);

% Determination of mode parameter
[OptimumParams, ~, P_std, C_std, T_std] = ecg_beat_fitter_25(meanECG,meanPhase);
% Ensure that inital bi estimates are positive
OptimumParams(8:14) = abs(OptimumParams(8:14));

% init state and its covariance
x0 = [phase(1); signal(1); signal(1); signal(1); OptimumParams];
P0 = 0.1*eye(25);

% State bound constraints
xmin = -inf * ones(size(x0));
xmax = inf * ones(size(x0));
% prevent bi from getting too small, which causes the linearized A matrix
% to explode and yields outliers during filtering
xmin(12:18) = 0.02;

% process and measurement noise
f = fs./diff(rWaveIndex);
omega_mean = 2*pi*mean(f);
omega_std = 2*pi*std(f);
w = [omega_mean, zeros(1,24)];
v = zeros(4,1);

% observation Y consisting of linear phase, p-wave, qrs-complex and t-wave
[P_signal, C_signal, T_signal] = ecg_windowing(signal, phase, -pi/6, pi/6, 30);
Y = [phase(:)'; P_signal(:)'; C_signal(:)'; T_signal(:)'];

[~, activity_only] = SNRmeasured(signal, rWaveIndex, fs);
emg_var = var(activity_only);

% Attempt at automatically tuning the noise covariance matrices based on
% some fundamental signal statistics.
% CAUTION, many manually tuned magic numbers ahead! Finding an optimal
% noise matrix tuning is _the_ crux of this method. If you find a better
% tuning method that works across many datasets, please do let me know.

% process noise covariance matrix
Q = eye(25);
Q(1,1) = omega_std * omega_std;
Q(2,2) = min(3e-3 * max(abs(P_signal)), 0.2*emg_var); % random noise on P signal
% This is typically way larger than EMG variance, hence it is already clear that this effectively just gates.
% Decreasing this variance leads to huge remaining QRS artifacts, because the waveform is not modeled sufficiently precisely by the Gaussian mixture model.
Q(3,3) = 3e-2 * max(abs(C_signal));  % random noise on QRS signal
Q(4,4) = min(3e-3 * max(abs(T_signal)), 0.2*emg_var); % random noise on T signal
Q(5:end,5:end) = diag([[2e-2; 2e-2; 6e-2; 6e-2; 6e-2; 2e-2; 2e-2] .* abs(OptimumParams(1:7)); % ai
	2e-2 * abs(OptimumParams(8:14));  % bi
	2e-3 * abs(OptimumParams(15:21))]); % thetai

% measurement noise covariance matrix
R = diag([omega_std^2, emg_var, emg_var, emg_var]);  

[Xhat, Xsmoothed, ~] = kalman_filter_25_sequential(Y, x0, P0, Q, R, w(:), v, fs, 0, xmin, xmax);

% Reconstruction of denoised ECG out of estimated p-wave, qrs-complex and t-wave
denoised_EKF = Xhat(2,:) + Xhat(3,:) + Xhat(4,:);
denoised_EKF = denoised_EKF(:);
denoised_EKS = Xsmoothed(2,:) + Xsmoothed(3,:) + Xsmoothed(4,:);
denoised_EKS = denoised_EKS(:);

EKF25 = signal - denoised_EKF;
EKS25 = signal - denoised_EKS;
end