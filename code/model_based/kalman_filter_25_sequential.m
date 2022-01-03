function [Xplus, Xsmoothed, gamma] = kalman_filter_25_sequential(Y, x0, P0, Q, R, w_mean, v_mean, fs, varargin)
% KALMAN_FILTER_25_SEQUENTIAL Implementation of sequential Extended Kalman
% Filter and Smoother used by EKF25. Uses sequential Modified Bryson-Frazier (MBF)
% smoother implementation without any matrix inversions during filtering or
% smoothing.
% If lower/upper bounds on states are provided, performs estimate
% projection to ensure xmin <= x <= xmax at all times.
%
%   INPUT:  Y           ->  rxn matrix where r is the number of measurement and
%                           n is the length of the measured signals. First row
%                           contains the phase measurement. Second to fourth
%                           row contains p-wave, qrs-complex and t-wave
%                           measurement, respectively.
%           x0          ->  Vector of length m with initial state.
%           P0          ->  mxm matrix with initial state covariance.
%           Q           ->  mxm matrix with process noise covariance.
%           R           ->  rxr matrix with oberservation noise covariance.
%           w_mean      ->  Vector of length m with process noise expection.
%           v_mean      ->  Vector of length r with oberservation noise
%                           expection.
%           fs          ->  Sampling frequency
%           use_mbf     ->  Choose algorithm: if 0, use sequential filtering and RTS. If 1 (default), use MBF.
%           xmin        ->  Lower bound for state values. (Optional)
%           xmax        ->  Upper bound for state values. (Optional)
%   OUTPUT: Xhat        ->  mxn matrix with estimated states by extended
%                           kalman filter
%           Xsmoothed   ->  mxn matrix with estimated states by extended
%                           kalman smoother.
%           gamma       ->  Innovation monitoring
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

if nargin < 9
    use_mbf = true;
else
    use_mbf = varargin{1};
end

if nargin < 10
    xmin = -inf * ones(size(x0));
else
    xmin = varargin{2};
end

if nargin < 11
    xmax = inf * ones(size(x0));
else
    xmax = varargin{3};
end

persistent dt;

dt = 1/fs;

n = length(Y);
m = length(x0);
r = size(Y, 1);

P_minus = P0;
x_minus = x0(:);

tic;

Xplus = zeros(m,n);
Xminus = zeros(m,n);
Pminus = zeros(m,m,n);

% Initialize innovation monitoring
gamma = zeros(4, n);
Ngamma = 50;
innov = zeros(Ngamma, r);
innov_covs = zeros(r, n);
innov_cov = zeros(r, r);

if use_mbf
    Vs = zeros(m, m, n);
    ydiffs = zeros(r, n);
else
    Pplus = zeros(m,m,n);
end

ydiff = zeros(r,1);
for i = 1:n
    
    Xminus(:,i) = x_minus;
    Pminus(:,:,i) = P_minus';
    
    %% UPDATE
    [M,N] = linearization_25(x_minus, w_mean, 1, dt); % linearisation of observation function
    assert(isdiag(N));  % required for the simple sequential implementation given here. 
    % For the more complex case see Dan Simon, optimal state estimation, p. 153.
    
    x_plus = x_minus;
    P_plus = P_minus;
    
    % Measurement noise covariance is diagonal! We can do sequential
    % updates.
    for j = 1:r

        innov_cov(j,j) = M(j,:) * P_plus * M(j,:)' + N(j,:) * R(j,j) * N(j,:)';

        if ~use_mbf
            % fully sequential update; does not play nicely with MBF
            Kj = P_plus * M(j,:)' / innov_cov(j,j);

            ydiff(j) = Y(j,i) - M(j,:) * x_plus + v_mean(j);
            if j == 1
                ydiff(1) = mod(ydiff(1) + pi, 2*pi) - pi;
            end

            x_plus = x_plus + Kj * ydiff(j);

            % Joseph stabilized Kalman covariance matrix
            P_plus = (eye(m) - Kj * M(j,:)) * P_plus * (eye(m) - Kj * M(j,:))' + Kj * N(j) * R(j,j) * N(j)' * Kj';
        end
    end
    
    if use_mbf
        K = P_plus * M' * diag(1./diag(innov_cov));  % inversion is cheap and stable since diagonal
        V = eye(m) - K * M;
        % Joseph stabilized Kalman covariance matrix
        P_plus = V * P_plus * V' + K * N * R * N' * K';
        Vs(:, :, i) = V;
        ydiff = Y(:,i) - M * x_plus + v_mean;
        ydiff(1) = mod(ydiff(1) + pi, 2*pi) - pi;
        ydiffs(:, i) = ydiff;
        x_plus = x_plus + K * ydiff;
    end
    
    % Enforce symmetry of covariance matrix
    P_plus = (P_plus + P_plus') / 2;
    
    % State projection to ensure keeping constraints satisfied.
    % This implements the approach detailed in section 2.4 of
    % "Kalman filtering with state constraints: a survey of linear and
    % nonlinear algorithms" by Dan Simon, IET Control Theory and
    % Applications, 2009.
    x_plus = min([xmax, max([xmin, x_plus], [], 2)], [], 2);
    
    Xplus(:,i) = x_plus;
    Pplus(:,:,i) = P_plus';
	
    %% PREDICTION
    x_minus = statePropagation(x_plus, w_mean, dt);
    
    [A,B] = linearization_25(x_plus, w_mean, 0, dt); % linearisation of process function
    
    P_minus = A * P_plus * A' + B * Q * B';
    
    % Monitor innovation signals
    innov(1:Ngamma-1, :) = innov(2:Ngamma, :);
    innov(Ngamma, :) = ydiff;
    innov_covs(:, i) = diag(innov_cov);  % we only need to store the diagonal
    for jj=1:4
        % (Time-smoothed) ratio of empirical variance of innovation signal
        % and expected variance. See "A Nonlinear Bayesian Filtering
        % Framework for ECG Denoising" (Sameni et al., 2007), eq. (13) for
        % details.
        gamma(jj, i) = ...
            sum(innov(max(1, Ngamma-i+1):end, jj).^2 ./ innov_covs(jj, max(1, i-Ngamma+1):i)') / Ngamma;
    end
end

Xsmoothed = zeros(m,n);
Psmoothed = zeros(m,m,n);

Psmoothed(:,:,n) = Pplus(:,:,n);
Xsmoothed(:,n) = Xplus(:,n);

if ~use_mbf
    % Rauch-Tung-Striebel-type smoothing, see e.g. S?rkk?
    % The required matrix inversion can lead to instabilities, hence this
    % is discouraged
    for i = n - 1:-1:1
        A = linearization_25(Xplus(:,i), w_mean, 0, dt);
        S = Pplus(:,:,i) * A' / (Pminus(:,:,i+1));
        Xsmoothed(:,i) = Xplus(:,i) + S * (Xsmoothed(:,i + 1) - Xminus(:,i + 1));
		
		% State projection to ensure keeping constraints satisfied.
        % This implements the approach detailed in section 2.4 of
        % "Kalman filtering with state constraints: a survey of linear and
        % nonlinear algorithms" by Dan Simon, IET Control Theory and
        % Applications, 2009.
        Xsmoothed(:,i) = min([xmax, max([xmin, Xsmoothed(:, i)], [], 2)], [], 2);
        
        Psmoothed(:,:,i) = Pplus(:,:,i) - S * (Pminus(:,:,i+1) - Psmoothed(:,:,i+1)) * S';
    end
else
    % Modified Bryson-Frazier Smoother (more stable than RTS: does not require inversion
    % of potentially near-singular P_minus)
    % Implemented after R. Gibbs (2011), "Square Root Modified Bryson?Frazier Smoother". 
    % (We do not use the square root version here.)
    Lambda_hat = zeros(size(A));
    lambda_hat = zeros(size(x0));
    for i = n:-1:1
        [A, ~] = linearization_25(Xplus(:, i), w_mean, 0, dt); % linearisation of process function
        [M, ~] = linearization_25(Xminus(:, i), w_mean, 1, dt); % linearisation of observation function

        % Smoothed covariance estimate
        Lambda_tilde = M' * diag(1./innov_covs(:, i)) * M + Vs(:, :, i)' * Lambda_hat * Vs(:, :, i);
        Lambda_hat = A' * Lambda_tilde * A;
        Psmoothed(:,:,i) = Pminus(:,:,i) - Pminus(:,:,i) * Lambda_tilde * Pminus(:,:,i);

        % Smoothed state estimate
        lambda_tilde = -M' * diag(1./innov_covs(:, i)) * ydiffs(:, i) + Vs(:, :, i)' * lambda_hat;
        lambda_hat = A' * lambda_tilde;
        Xsmoothed(:,i) = Xminus(:, i) - Pminus(:,:,i) * lambda_tilde;
		
        % State projection to ensure keeping constraints satisfied.
        % This implements the approach detailed in section 2.4 of
        % "Kalman filtering with state constraints: a survey of linear and
        % nonlinear algorithms" by Dan Simon, IET Control Theory and
        % Applications, 2009.
        Xsmoothed(:,i) = min([xmax, max([xmin, Xsmoothed(:, i)], [], 2)], [], 2);
    end
end

end


function x_minus = statePropagation(x_plus, w_mean, dt)

phi = x_plus(1);
P = x_plus(2);
C = x_plus(3);
T = x_plus(4);
alphai = x_plus(5:11);
bi = x_plus(12:18);
thetai = x_plus(19:25);

omega = w_mean(1);

x_minus = zeros(size(x_plus));

% phi_k+1 = (phi_k + omega * dt) mod (2 * pi)
x_minus(1) = phi + omega * dt;                                                    
x_minus(1) = mod(x_minus(1) + pi, 2*pi) - pi;

%dtetai = mod(x_minus(1) - thetai, pi);
dtetai = mod(x_minus(1) - thetai + pi, 2*pi) - pi;

% P_k+1
x_minus(2) = -dt * sum(alphai(1:2) * omega ./ (bi(1:2) .^2) .* dtetai(1:2) ...
    .* exp(-dtetai(1:2) .^2 ./ (2 * bi(1:2) .^2))) + P + w_mean(2);
% C_k+1
x_minus(3) = -dt * sum(alphai(3:5) * omega ./ (bi(3:5) .^2) .* dtetai(3:5) ...
    .* exp(-dtetai(3:5) .^2 ./ (2 * bi(3:5) .^2))) + C + w_mean(3);
% T_k+1
x_minus(4) = -dt * sum(alphai(6:7) * omega ./ (bi(6:7) .^2) .* dtetai(6:7) ...
    .* exp(-dtetai(6:7) .^2 ./ (2 * bi(6:7) .^2))) + T + w_mean(4);

% alphai_k+1, bi_k+1, thetai_k+1
x_minus(5:end) = x_plus(5:end) + w_mean(5:end);

end