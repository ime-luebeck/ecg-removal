function [P1, f] = calcFrequencySpectrum(signal, varargin)
%CALCFREQUENCYSPECTRUM 
%
%   INPUT:  signal  ->  Vector with one signal.
%           fs      ->  Sampling frequency in Hz.
%   OUPUT:  P1      ->  Vector with spectral power density estimate of 'signal'
%                       calculated by multi-taper method
%           f       ->  Vector with same length as 'P1' that contains
%                       corresponding frequency scale.
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

if nargin <= 1
    fs = 1024;
else 
    fs = varargin{1};
end

% Use the default Hamming window, 50% overlap, and default choice of NFFT
%[P1, f] = pwelch(signal, window_len, [], nfft, fs);
[P1, f] = pmtm(signal, 15, 4096, fs);

end