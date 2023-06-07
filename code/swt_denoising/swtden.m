function [out, wDEC, thresholds, gates] = swtden(x, peaks, fs, varargin)
% SWTDEN Shrinkage Denoising using a-trous wavelet decomposition (swt) 
%
% Copyright 2021 Institute for Electrical Engineering in Medicine, 
% University of Luebeck
% Jan Graßhoff
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

x = x(:)';
peaks = peaks(:)';

% Use soft ('s') or hard ('h') thresholding
if nargin > 3
	sorh = varargin{1};
else
	sorh = 'h';
end

% Number of downsampling steps
if nargin > 4
	n = varargin{2};
else
	n = max(2, round(log2(fs)-6));  % This yields n=3 for fs=512, n=4 for fs=1024, n=5 for fs=2048, etc.
end

% Wavelet to use
if nargin > 5
	w = varargin{3};
else
	w = 'db4';
end

% thresholds: above how many standard deviations from the mean is a signal
% no longer considered 'noise'?
if nargin > 6
	fixed_thresh = varargin{4};
else
	fixed_thresh = 4.5;
end

% Calculate gates
gates = get_gates(peaks, round(0.1*fs));

% Signal Extension by zero padding
pow = 2^n;
l = length(x);
len_extended = ceil(l/pow) * pow;
zpd = zeros(len_extended-l, 1);
x_zpd = [x; zpd];
gates = [gates(:); zpd];

% Wavelet decomposition of x.
wDEC = swt(x_zpd, n, w);

% Gate out R-peaks in wavelet subbands
wDEC_gated = wDEC;
wDEC_gated(:, gates==1) = NaN;

% Custom threshold coefficients
win_len = fs;
s = noise_est(wDEC_gated(1:end-1, :), win_len);

thresholds = zeros(size(wDEC));
wxd = wDEC;

for k = 1:(n)
    threshold = fixed_thresh * s(k, :);
    thresholds(k, :) = threshold;
    wxd(k, :) = wthresh(wDEC(k, :), sorh, threshold);
end

% Wavelet reconstruction
xd = iswt(wxd, w);

% Return results
wDEC = wDEC(:, 1:l);
xd = xd(1:l);
thresholds = thresholds(:, 1:l);
gates = gates(1:l)';

out = x - xd;

end


function stdc = noise_est(signal, win_len)
% Estimate noise level

nblev = size(signal, 1);
stdc = zeros(size(signal));

for k = 1:nblev
    % robust (moving) standard deviation estimator
    stdc(k, :) = medfilt1(abs(signal(k, :)), win_len, 'omitnan') / 0.6745;
    stdc(k, 1:ceil(win_len/2)) = stdc(k, ceil(win_len / 2));
    stdc(k, end - ceil(win_len/2):end) = stdc(k, end - ceil(win_len/2));
end

end


function gates = get_gates(rpeaks, gate_width_samples)

gate_width_samples = floor(gate_width_samples/2) * 2;

rpeak_indx = find(rpeaks==1);

gates = zeros(size(rpeaks));
for i=1:length(rpeak_indx)
    gates(max(rpeak_indx(i) - gate_width_samples/2, 1):min(rpeak_indx(i) + gate_width_samples/2, length(rpeaks))) = 1;
end

end

