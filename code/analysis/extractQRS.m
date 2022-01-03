function [QRScomplexes, QRS_mask] = extractQRS(signal, rWaveIndex, fs)
%EXTRACTQRS Sectioning of input signal with the help of detected r peaks.
%   QRS-complexes and part of t-Wave are selected from an ECG signal or
%   cardiac contaminated EMG signal.
%   INPUT:  signal      ->  Vector with EMG or ECG signal.
%           rWaveIndex  ->  Vector that contains indices of detected r
%                           peaks in 'signal'.
%           fs          ->  Sampling frequency
%   OUTPUT: QRScomplexes->  Vector that contains only selected QRS
%                           complexes
%           QRS_mask    ->  Vector with same length as 'signal' that is 1 if
%                           sample belongs to QRS complexes, 0 else.
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

numQRS = length(rWaveIndex);
% Mean QRS complex duration is 80 ms
meanQRSlengthLeft = ceil(0.15 * fs);
meanQRSlengthRight = ceil(0.15 * fs);
startQRS = max(rWaveIndex-meanQRSlengthLeft, 1);
endQRS = min(rWaveIndex+meanQRSlengthRight, length(signal));

QRS_mask = zeros(size(signal));
for i = 1:numQRS
    QRS_mask(startQRS(i):endQRS(i)) = 1;
end

QRScomplexes = signal(QRS_mask~=0);

end