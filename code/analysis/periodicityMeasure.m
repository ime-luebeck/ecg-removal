function pm = periodicityMeasure(signal, rWaveIndex, phase)
%PERIODICITYMEASURE Calculates the ECG-related periodicity measure
%proposed by Sameni, Shamsollahi and Jutten (2008).
%   Calculates Pearson's correlation coefficient between a signal
%   and a version of the signal that is shifted adaptively, such that
%   each sample is compared to the sample in the following ECG beat
%   that occurs at the same relative phase of the ECG cycle.
%
%   Returns a number between 0 (no correlation) and 1 (exact
%   correlation) that indicates the relative strength of the ECG 
%   component in the signal.
%
%   Has been proposed in Sameni, Shamsollahi, Jutten (2008): "Model-
%   based Bayesian filtering of cardiac contaminants from biomedical
%   recordings."
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

m = length(rWaveIndex);

kplustau = zeros(size(signal));

start = 40;
currentPhase = phase(rWaveIndex(start):rWaveIndex(start+1)-1);

for i = start+1:m-1 % go through phase sections
    nextPhase = phase(rWaveIndex(i):rWaveIndex(i+1)-1);
    currentTau = findTau(currentPhase, nextPhase);
    kplustau(rWaveIndex(i-1):rWaveIndex(i)-1) = rWaveIndex(i) + currentTau;
    currentPhase = nextPhase;
end

kplustau = kplustau(rWaveIndex(start):rWaveIndex(m-1)-1);
k = rWaveIndex(start):rWaveIndex(m-1)-1;

corrcoeffs = abs(corrcoef(signal(k), signal(kplustau)));
pm = corrcoeffs(1,2);

end


function tau = findTau(currentPhase, nextPhase)
%FINDTAU Used by function periodicityMeasure().
%   INPUT:  currentPhase    ->  Vector that contains one phase section
%                               between two r peaks (values between -pi and
%                               pi).
%           nextPhase       ->  Vector with following phase section. Both
%                               vectors don't need to have same length.
%   OUTPUT: tau             ->  Vector with same length as 'currentPhase'.
%                               tau defines the index of dual sample in 'nextphase'
%                               for every sample of 'currentphase'.

tau = zeros(size(currentPhase));

for i = 1:length(currentPhase)
    [~, tau(i)] = min(abs(nextPhase - currentPhase(i)));
end

end
