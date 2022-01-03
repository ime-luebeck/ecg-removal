function [P_signal, C_signal, T_signal] = ecg_windowing(ECGsignal, phase, varargin)
%ECGWINDOWING Windowing of the ECG signal to seperate into three states: P-wave, QRS-complex and T-wave
%   The windows are formed by sigmoidal functions. 
%   INPUT:  ECGsignal
%           phase: corresponding (linear or non-linear) phase
%           a: boundary between P-wave and QRS-complex, -pi/6 by default
%           b: boundary between QRS-complex and T-wave, pi/6 by default
%           gamma: shape of the sigmoidal function, 30 by default
%   OUTPUT: P_signal: observed P-state
%           C_signal: observed QRS-state
%           T_signal: observed T-state
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

if nargin == 2
    a = -pi/6;
    b = pi/6;
    gamma = 30;
elseif nargin == 5
    a = varargin{1};
    b = varargin{2};
    gamma = varargin{3};
end

Pwindow = max(0, sigmf(mod(phase+pi+pi, 2*pi) - pi, [gamma, 0]) ...
               - sigmf(mod(phase+pi-a, 2*pi) - pi, [gamma, 0]));
Cwindow = max(0, sigmf(mod(phase+pi-a, 2*pi) - pi, [gamma, 0]) ...
               - sigmf(mod(phase+pi-b, 2*pi) - pi, [gamma, 0]));
Twindow = max(0, sigmf(mod(phase+pi-b, 2*pi) - pi, [gamma, 0]) ...
               - sigmf(mod(phase+pi-pi, 2*pi) - pi, [gamma, 0]));

% Windows should sum to 1!
assert(all(abs(Pwindow+Cwindow+Twindow - 1) < 1e-5))
           
P_signal = ECGsignal(:).*Pwindow(:);
C_signal = ECGsignal(:).*Cwindow(:);
T_signal = ECGsignal(:).*Twindow(:);

end

