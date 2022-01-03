function [X_filt, X_smooth, P_filt, P_smooth, MeasNoiseFilt, MeasNoiseSmoothed, innovation, energy, ...
    dEnergydQroot, dEnergydRroot] = KF_RTS(Y, A, C, Q, R, varargin)
% KF_RTS sequential Kalman Filter and Rauch-Tung-Striebel Smoother
% Linear Kalman Filter and Smoother; sequential implementation of Joseph
% stabilized, symmetrized form. Only works for diagonal measurement noise
% covariance matrices R!
%
% To use a time-varying system description, pass 3d arrays for A, C, and/or Q,
% where the last index iterates over sample indices and/or an array for R,
% where the rows indicate different samples.
%   x(k+1) = A(k) * x(k) + q  [ + B * u(k) ]
%   y(k)   = C(k) * x(k) + r
%
% For details on the KF and its sequential implementation and the RTS see Dan Simon, Optimal State Estimation (2006), p. 150.
% The formulation of the energy function calculation is due to Simo Särkkä, "Bayesian filtering and smoothing".
% The formulation of state constraints is due to section 2.4 of "Kalman filtering with state constraints: a survey of linear and
% nonlinear algorithms" by Dan Simon, IET Control Theory and Applications, 2009.
% For an explanation of what the weighted version of this does (and why you might want it), see
% Eike Petersen, "Model-based Probabilistic Inference for Monitoring Respiratory Effort using the Surface Electromyogram"
% (Dissertation, forthcoming).
% Some of the implementation details are my own. (Eike Petersen)
%
%   INPUT:  Y           ->  rxN matrix where r is the number of measurement and
%                           N is the length of the measured signals.
%           A           ->  System matrix.
%           C           ->  Observation matrix
%           Q           ->  Process noise covariance.
%                           Should be either of the following:
%                           - mxm, specifying a constant matrix.
%                           - mxmxN, specifying a time-varying matrix.
%                           - 1xm, specifying the diagonal of a constant
%                           matrix.
%                           - Nxm, specifying a time-varying diagonal
%                           matrix.
%           R           ->  1xr vector with oberservation noise covariance
%                           matrix diagonal entries. If time-varying, this
%                           should be Nxr.
%           x0          ->  Vector of length m with initial state. (Optional, positional)
%           P0          ->  mxm matrix with initial state covariance. (Optional, positional)
%           handle_nans ->  handle 'nan' values in measurements or
%                           observation matrices by simply treating this
%                           measurement as completely unknown (infinite 
%                           measurement variance). Default is true. (Optional, name-value)
%           xmin        ->  Lower bound for state values. (Optional, name-value)
%           xmax        ->  Upper bound for state values. (Optional, name-value)
%           sample_weights -> Weight factors for individual measurement samples.
%                           1xN matrix. The default are equal weights=1.
%                           If provided, biases the estimation towards samples
%                           with higher weight but also increases estimator 
%                           variance. (Optional, name-value)
%           Qroot       ->  Matrix square root of the (in this case, necessarily static) process
%                           noise covariance matrix Q=Qroot'*Qroot. Required for calculating
%                           the derivative of the average log likelihood wrt exactly this matrix. (Optional, name-value)
%           Rroot       ->  Matrix square root of the (in this case, necessarily static) measurement
%                           noise covariance matrix R=Rroot'*Rroot. Required for calculating
%                           the derivative of the average log likelihood wrt exactly this matrix. (Optional, name-value)
%           B           ->  m x o, where o is the number of input signals. (Optional, name-value)
%           u           ->  N x o, input signals u(k). (Optional, name-value)
%           Qmean       ->  Process noise mean vector. If not provided, Qmean=0 is assumed. (Optional, name-value)

%
%   OUTPUT: X_filt            ->  mxN matrix with estimated states by
%                                 Kalman filter
%           X_smooth          ->  mxN matrix with estimated states by 
%                                 Rauch-Tung-Striebel smoother.
%           P_filt            ->  Covariance matrices of filtered state signal (mxmxN)
%           P_smooth          ->  Covariance matrices of smoothed state signal (mxmxN)
%           MeasNoiseFilt     ->  Filtered measurement error estimates
%           MeasNoiseSmoothed ->  Smoothed measurement error estimates
%           innovation        ->  Time course of the innovation signal / filter prediction errors
%           Energy            ->  energy = -log p(y)
%                                 This is a measure of model plausibility and can be used to compare different models
%                                 or data sets.Samples with NaN-values or R=inf are ignored for this calculation.
%           dEnergydQroot       ->  Derivative of the energy function described above wrt Qroot,
%                                 i.e., the square root of the process noise covariance matrix.
%           dEnergydRroot       ->  Derivative of the energy function described above wrt Rroot
%                                 i.e., the square root of the measurement noise covariance matrix.
%
%
%   Eike Petersen, 2015-2021

%% INPUT HANDLING

p = inputParser;

% Positional arguments
addRequired(p, 'Y', @ismatrix);
addRequired(p, 'A', @(A) isnumeric(A) && ndims(A) <= 3 && ~any(isnan(A(:))));
addRequired(p, 'C', @(C) isnumeric(C) && ndims(C) <= 3 && ~all(isnan(C(:))));
addRequired(p, 'Q', @(Q) isnumeric(Q) && ndims(Q) <= 3 && ~any(isnan(Q(:))));
addRequired(p, 'R', @(R) isnumeric(R) && ndims(R) <= 3);
addOptional(p, 'x0', [], @(x0) isnumeric(x0) && ismatrix(x0) && ~any(isnan(x0)));
addOptional(p, 'P0', [], @(P0) isnumeric(P0) && ismatrix(P0) && ~any(isnan(P0(:))));
% Name-value arguments
addParameter(p, 'handle_nans', true, @islogical);
addParameter(p, 'xmin', [], @isnumeric);
addParameter(p, 'xmax', [], @isnumeric);
addParameter(p, 'sample_weights', [], @(w) isnumeric(w) && ~any(isnan(w(:))));
addParameter(p, 'Qroot', [], @(Qr) isnumeric(Qr) && ~any(isnan(Qr(:))));
addParameter(p, 'Rroot', [], @(Rr) isnumeric(Rr) && ~any(isnan(Rr(:))));
addParameter(p, 'B', [], @(B) isnumeric(B) && ismatrix(B) && ~any(isnan(B(:))));
addParameter(p, 'u', [], @(u) isnumeric(u) && ismatrix(u) && ~any(isnan(u(:))));
addParameter(p, 'Qmean', [], @(Qm) isnumeric(Qm) && ismatrix(Qm) && ~any(isnan(Qm(:))));

parse(p, Y, A, C, Q, R, varargin{:})
handle_nans = p.Results.handle_nans;

if isempty(p.Results.x0)
    x0 = zeros(size(A, 1), 1);
else
    x0 = p.Results.x0;
    assert(isvector(x0));
end

if isempty(p.Results.P0)
    P0 = 1e10 * eye(size(A, 1));
else
    P0 = p.Results.P0;
    assert(ismatrix(P0));
end

if isempty(p.Results.xmin)
    xmin = -inf * ones(size(x0));
else
    xmin = make_col_vec(p.Results.xmin);
end

if isempty(p.Results.xmax)
    xmax = inf * ones(size(x0));
else
    xmax = make_col_vec(p.Results.xmax);
end

m = length(x0);
[r, N] = size(Y);

B = p.Results.B;
u = p.Results.u;

if isempty(p.Results.Qmean)
    Qmean = zeros(m, 1);
else
    Qmean = p.Results.Qmean;
end

Qroot = p.Results.Qroot;
Rroot = p.Results.Rroot;

if isempty(p.Results.sample_weights)
    sample_weights = ones(N, 1);
    non_standard_weights = false;
else
    sample_weights = p.Results.sample_weights;
    if all(sample_weights == 1)
        non_standard_weights = false;
    else
        non_standard_weights = true;
        if nargout >= 9
            warning('Energy derivatives are only implemented for standard weights and constant matrices, will return NaN.');
        end
    end
end
assert(abs(sum(sample_weights) - N) < 1)


%% INITIALIZATION

% State transition matrix
function Ai = get_Ai(ii)
    if ndims(A) == 3
        Ai = A(:, :, ii);
    else
        Ai = A;
    end
end

% Measurement matrix
function Ci = get_Ci(ii)
    if ndims(C) == 3
        Ci = C(:, :, ii);
        assert(~any(isinf(Ci(:))))
    else
        Ci = C;
        assert(~any(isinf(Ci(:))))
    end
end

% Measurement noise covariance
assert(size(R, 2) == r);
assert(size(R, 1) == 1 || size(R, 1) == N);

function Ri = get_Ri(ii)
    if size(R, 1) > 1
        Ri = diag(R(ii, :));
    else
        Ri = diag(R);
    end
    % For any unknown entries, assume the worst
    Ri(isnan(Ri)) = inf;
end

% Process noise covariance
if ndims(Q) == 1
    % single diagonal
    assert(length(Q) == m);
elseif ismatrix(Q)
    assert(size(Q, 1) == 1 || size(Q, 1) == m || size(Q, 1) == N);
    assert(size(Q, 2) == m);
else
    assert(ndims(Q) == 3)
    assert(size(Q, 3) == N);
    assert(size(Q, 1) == m && size(Q, 2) == m);
end

function Qi = get_Qi(ii)
    if ndims(Q) == 3
        Qi = Q(:, :, ii);
    else
        if all(size(Q) == [m, m])
            Qi = Q;
        elseif all(size(Q) == [1, m])
            Qi = diag(Q);
        elseif all(size(Q) == [N, m])
            Qi = diag(Q(ii, :));
        end
    end
end

if ~isempty(Qroot)
    % verify that the user-provided matrix square root is actually a reasonable square root
    QQ = Qroot'*Qroot;
    Q1 = get_Qi(1);
    assert(all(abs((QQ(:) - Q1(:))) < 1e-6));
end
if ~isempty(Rroot)
    % verify that the user-provided matrix square root is actually a reasonable square root
    RR = Rroot'*Rroot;
    R1 = get_Ri(1);
    assert(all(abs((RR(:) - R1(:))) < 1e-6));
end

P_minus_ii = P0;
P_minus_ii_w = P0;
x_minus = x0(:);

X_filt = zeros(m, N);
P_plus = zeros(m, m, N);
X_minus = zeros(m, N);
P_minus = zeros(m, m, N);

if non_standard_weights
    P_minus_w = zeros(m, m, N);
end

if nargout > 2
    P_filt = zeros(m, m, N);
end

if nargout > 4
    MeasNoiseFilt = zeros(r, N);
end

innov_cov = zeros(r,r);
if non_standard_weights && nargout > 7
    innov_cov_w = zeros(r,r);
end

ydiff = zeros(r,1);

if nargout > 6
    innovation = zeros(r, N);
end

if nargout > 7
    % energy function from Särkkä book, phi(params) = -log p(y | params) - log p(params)
    % initialize with phi = -log p(params) = 0, i.e. assume a uniform prior.
    energy = 0;
    nan_or_inf_energy_steps = 0;
end

if nargout > 8
    if ~(max(size(Q)) == m)
        warning('Energy derivate with respect to Q is only implemented for constant Q. Returning NaN.')
        dEnergydQroot = nan;
        calc_dEnergydQroot = false;
    elseif non_standard_weights
        dEnergydQroot = nan;
        calc_dEnergydQroot = false;
    else
        dEnergydQroot = 0;
        if isempty(Qroot)
            % calculate square root of Q
            Qroot = sqrtm(Q);
        end
        Q1 = get_Qi(1);
        calc_dEnergydQroot = true;
    end
else
    calc_dEnergydQroot = false;
end

if nargout > 9
    if ~ismatrix(C)
        warning('Can only calculate the derivative of the energy function w.r.t. measurement noise covariance if C is constant. Returning NaN.')
        dEnergydRroot = nan;
        calc_dEnergydRroot = false;
    elseif non_standard_weights
        dEnergydRroot = nan;
        calc_dEnergydRroot = false;
    else
        RTS_D = 0;
        calc_dEnergydRroot = true;
    end
else
    calc_dEnergydRroot = false;
end


%% FILTERING

for ii = 1:N

    Ai = get_Ai(ii);
    Ci = get_Ci(ii);
    Qi = get_Qi(ii);
    Ri = get_Ri(ii);
    
    X_minus(:,ii) = x_minus;
    P_minus(:,:,ii) = P_minus_ii';
    x_plus = x_minus;
    P_plus_ii = P_minus_ii;

    if non_standard_weights
        P_plus_ii_w = P_minus_ii_w;
        P_minus_w(:,:,ii) = P_minus_ii_w';
    end
    
    % -----------
    % UPDATE STEP
    % -----------
    
    if nargout > 6
        innovation(:, ii) = Y(:, ii) - Ci * x_plus;
    end
    
    % Sequential updates
    for j = 1:r
        
        innov_cov(j,j) = Ci(j,:) * P_plus_ii * Ci(j,:)' + Ri(j,j) / sample_weights(ii);

        if non_standard_weights && nargout > 7
            innov_cov_w(j,j) = Ci(j,:) * P_plus_ii_w * Ci(j,:)' + Ri(j,j);
        end        
        
        K = P_plus_ii * Ci(j,:)' / innov_cov(j,j);
        ydiff(j) = Y(j,ii) - Ci(j,:) * x_plus;
        
        % Output NaN handling
        if handle_nans && (any(isnan(K)) || isnan(ydiff(j)))
            % Something is unknown, hence assume Ri -> inf => K -> 0.
            K = zeros(size(K));
            ydiff(j) = 0; % value doesn't matter
            KC = zeros(size(K*Ci(j, :))); % assuming Ci(j, :) < inf; see assertion above
        else
            KC = K * Ci(j, :);
        end
        
        x_plus = x_plus + K * ydiff(j);
        
        % Joseph stabilized Kalman covariance matrix
        if ~isinf(Ri(j, j) / sample_weights(ii))
            P_plus_ii = (eye(m) - KC) * P_plus_ii * (eye(m) - KC)' + K * Ri(j,j) * K' / sample_weights(ii);
        else
            % The last summand K*Ri(j,j)*K' is what causes problems: this is 0*inf*0=nan.
            % However, the analytical limit of this term is 0, so we just omit it.
            P_plus_ii = (eye(m) - KC) * P_plus_ii * (eye(m) - KC)';
        end

        % Enforce symmetry of covariance matrix
        P_plus_ii = (P_plus_ii + P_plus_ii') / 2;
        
        if non_standard_weights
            % Joseph stabilized Kalman covariance matrix
            if ~isinf(Ri(j, j))
                P_plus_ii_w = (eye(m) - KC) * P_plus_ii_w * (eye(m) - KC)' + K * Ri(j,j) * K';
            else
                % The last summand K*Ri(j,j)*K' is what causes problems: this is 0*inf*0=nan.
                % However, the analytical limit of this term is 0, so we just omit it.
                P_plus_ii_w = (eye(m) - KC) * P_plus_ii_w * (eye(m) - KC)';
            end
            P_plus_ii_w = (P_plus_ii_w + P_plus_ii_w') / 2;
        end
    end
    
    % State projection to ensure keeping constraints satisfied.
    % This implements the approach detailed in section 2.4 of
    % "Kalman filtering with state constraints: a survey of linear and
    % nonlinear algorithms" by Dan Simon, IET Control Theory and
    % Applications, 2009.
    x_plus = min([xmax, max([xmin, x_plus], [], 2)], [], 2);
    
    X_filt(:,ii) = x_plus;
    P_plus(:,:,ii) = P_plus_ii';
    
    if non_standard_weights
        P_filt(:, :, ii) = P_plus_ii_w';
    else
        P_filt(:, :, ii) = P_plus_ii';
    end
    
    if nargout > 4
        MeasNoiseFilt(:, ii) = Y(:, ii) - Ci * x_plus;
    end
    
    % ---------------
    % PREDICTION STEP
    % ---------------
    
    if isempty(B)
        x_minus = Ai * x_plus + Qmean;
        P_minus_ii = Ai * P_plus_ii * Ai' + Qi;
    else
        x_minus = Ai * x_plus + B * u(ii, :)' + Qmean;
        P_minus_ii = Ai * P_plus_ii * Ai' + Qi;
    end
    
    if non_standard_weights
        P_minus_ii_w = Ai * P_plus_ii_w * Ai' + Qi;
    end
    
    if nargout > 7
        if ~any(isnan(innovation(:, ii))) && ~any(isnan(innov_cov(:)))
            % update energy function: phi_k(params) = phi_{k-1}(params) - log p(y(k) | y(1:k-1), params)
            % [Eqs. 12.11+12.38 in Särkkä book]
            if ~non_standard_weights
                energy_ii = 0.5 * log(det(2*pi*innov_cov)) + 0.5 * innovation(:, ii)' * (innov_cov \ innovation(:, ii));
            else
                energy_ii = 0.5 * log(det(2*pi*innov_cov_w)) + 0.5 * innovation(:, ii)' * (innov_cov_w \ innovation(:, ii));
            end
            if ~(isnan(energy_ii) || isinf(energy_ii))
                energy = energy + sample_weights(ii) * energy_ii;
            else
                nan_or_inf_energy_steps = nan_or_inf_energy_steps + 1;
            end
        else
            nan_or_inf_energy_steps = nan_or_inf_energy_steps + 1;
        end
    end
    
    if calc_dEnergydRroot
        RTS_D = RTS_D + 1/N * (Y(:,ii) * Y(:,ii)');
    end
end

if nargout > 1
    %% SMOOTHING
    % ----------
    % Rauch-Tung-Striebel-type smoothing, see, e.g., Simo Särkkä (2013): "Bayesian Filtering and Smoothing"

    X_smooth = zeros(m, N);
    P_smooth_internal = zeros(m, m, N);
    
    % initialize last time step with filter results
    P_smooth_internal(:, :, N) = P_plus(:, :, N);
    X_smooth(:, N) = X_filt(:, N);

    if nargout > 3
        P_smooth = zeros(m, m, N);
        P_smooth(:, :, N) = P_filt(:, :, N); 
    end
    
    if nargout > 5
        MeasNoiseSmoothed = zeros(size(Y, 1), N);
        MeasNoiseSmoothed(:, N) = MeasNoiseFilt(:, N);
    end
    
    if calc_dEnergydQroot
        % initialize quantities required to evaluate the derivative of the energy function as defined by Särkkä
        RTS_C = 0;
        RTS_Sigma = 1/N * (P_smooth_internal(:, :, N) + X_smooth(:, N) * X_smooth(:, N)');
    end
    
    if calc_dEnergydRroot
        % initialize quantities required to evaluate the derivative of the energy function as defined by Särkkä
        RTS_B = 1/N * (Y(:, N) * X_smooth(:, N)');
    end
    
    % Set up warning handling to prevent spamming the command line in case of ill-conditioning
    warning('off', 'MATLAB:illConditionedMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    singular_matrix_warnings_found = 0;
    lastwarn('');
    for ii = N - 1:-1:1
        Ai = get_Ai(ii);
        Qi = get_Qi(ii);

        if ~any(Qi(:))
            % Qi == 0. In this case, P_minus_ii = P_{k+1}^- = A * P_{k-1}^T * A^T
            % and we have (analytically) that S = A^-1.
            S0 = inv(Ai);
        else
            % use standard text book formula
            S0 = P_plus(:,:, ii) * Ai' / P_minus(:,:, ii+1);
        end

        X_smooth(:, ii) = X_filt(:, ii) + S0 * (X_smooth(:, ii + 1) - X_minus(:, ii + 1));
        % State projection to ensure keeping constraints satisfied.
        % This implements the approach detailed in section 2.4 of
        % "Kalman filtering with state constraints: a survey of linear and
        % nonlinear algorithms" by Dan Simon, IET Control Theory and
        % Applications, 2009.
        X_smooth(:, ii) = min([xmax, max([xmin, X_smooth(:, ii)], [], 2)], [], 2);   

        P_smooth_internal(:, :, ii) = P_plus(:, :, ii) + S0 * (P_smooth_internal(:, :, ii+1) - P_minus(:, :, ii+1)) * S0';

        if nargout > 3
            if non_standard_weights
                P_smooth(:, :, ii) = P_filt(:, :, ii) + S0 * (P_smooth(:, :, ii+1) - P_minus_w(:, :, ii+1)) * S0';
            else
                P_smooth(:, :, ii) = P_smooth_internal(:, :, ii);
            end
        end
        
        if nargout > 5
            Ci = get_Ci(ii);
            MeasNoiseSmoothed(:, ii) = Y(:, ii) - Ci * X_smooth(:, ii);
        end

        if calc_dEnergydQroot
            % compute quantities required to evaluate the derivative of the energy function as defined by Särkkä
            RTS_C = RTS_C + 1/N * (P_smooth_internal(:, :, ii+1) * S0' + X_smooth(:, ii+1) * X_smooth(:, ii)');
            RTS_Sigma = RTS_Sigma + 1/N * (P_smooth_internal(:, :, ii) + X_smooth(:, ii) * X_smooth(:, ii)');
        end
        
        if calc_dEnergydRroot
            % compute quantities required to evaluate the derivative of the energy function as defined by Särkkä
            RTS_B = RTS_B + 1/N * (Y(:, ii) * X_smooth(:, ii)');
        end
        
        [~, msgidlast] = lastwarn;
        if strcmp(msgidlast,'MATLAB:illConditionedMatrix') || strcmp(msgidlast,'MATLAB:nearlySingularMatrix') ...
                || strcmp(msgidlast, 'MATLAB:singularMatrix')
            singular_matrix_warnings_found = singular_matrix_warnings_found + 1;
        else
            disp(msgidlast)
        end

        lastwarn('');
    end  % end smoother loop

    if singular_matrix_warnings_found > 0
        warning('kalmanFilterSeq:Smoothing:SingularMatrices', ...
            'Matrix close to singular or singular during %d smoothing steps.', singular_matrix_warnings_found);
    end
    warning('on', 'MATLAB:illConditionedMatrix');
    warning('on', 'MATLAB:nearlySingularMatrix');
    warning('on', 'MATLAB:singularMatrix');
    
    if calc_dEnergydQroot
        % Phi from Eq. 12.44 in Särkkä's book is identical to Sigma except for first / last summands
        RTS_Phi = RTS_Sigma - 1/N * (P_smooth_internal(:, :, N) + X_smooth(:, N) * X_smooth(:, N)');
        % In Särkkä's formulation, we should now add the smoothed estimate of the initial state
        % However, I don't have that here, because in my implementation there is a measurement
        % directly for the first time step. Thus we need to recover an artificial first state...
        Xinitial = Ai \ X_smooth(:, 1);
        Pinitial = Ai \ (P_smooth_internal(:, :, 1) + Q1) / Ai';
        RTS_Phi = RTS_Phi + 1/N * (Pinitial + Xinitial * Xinitial');
        P_plus_0 = Ai \ (P_minus(:, :, 1) - Qi) / Ai';
        S0 =  P_plus_0 * Ai' / P_minus(:, :, 1);
        RTS_C = RTS_C + 1/N * (P_smooth_internal(:, :, 1) * S0' + X_smooth(:, 1) * Xinitial');
        % derivative of the energy function wrt Q [Eqs A.17+A.18 in Särkkä book]
        dEnergydQroot = zeros(size(Qroot));
        % there's probably a more efficient vectorized / matrix formulation for this...
        for ii = 1:size(Qroot, 1)
            for jj = 1:size(Qroot, 1)
                % derivative of the proc noise covariance Q=L'L wrt entry (i,j) of L=Qroot
                Jij = zeros(size(Qroot));
                Jij(ii, jj) = 1;
                Jji = zeros(size(Qroot));
                Jji(jj, ii) = 1;                
                dQdQrootij = Qroot' * Jij + Jji * Qroot;
                % derivative of the energy function wrt L=Qroot
                dEnergydQroot(ii, jj) = ...
                    +0.5 * N * trace(Q1 \ dQdQrootij) ...
                    -0.5 * N * trace((Q1 \ dQdQrootij) * ...
                        (Q1 \ (RTS_Sigma - RTS_C * Ai' - Ai * RTS_C' + Ai * RTS_Phi * Ai')));
            end
        end
    end
    
    if calc_dEnergydRroot
        dEnergydRroot = zeros(size(Rroot));
        % there's probably a more efficient vectorized / matrix formulation for this...
        for ii = 1:size(Rroot, 1)
            for jj = 1:size(Rroot, 1)
                % derivative of the proc noise covariance Q=L'L wrt entry (i,j) of L=Qroot
                Jij = zeros(size(Rroot));
                Jij(ii, jj) = 1;
                Jji = zeros(size(Qroot));
                Jji(jj, ii) = 1;                    
                dRdRrootij = Rroot' * Jij + Jji * Rroot;
                % derivative of the energy function wrt L=Qroot
                dEnergydRroot(ii, jj) = ...
                    +0.5 * N * trace(R1 \ dRdRrootij) ...
                    -0.5 * N * trace((R1 \ dRdRrootij) * ...
                        (R1 \ (RTS_D - RTS_B * C' - C * RTS_B' + C * RTS_Sigma * C')));
            end
        end
    end        
end

end