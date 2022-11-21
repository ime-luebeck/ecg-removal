function [SNRimprovement, SNRinput] = analyzeSNRimprovement(dataWithRPeaks, fs, ...
    data, varargin)
%ANALYZESNRIMPROVMENT
%   INPUT:  dataWithRPeaks  ->  nx5 cell-matrix, where n is the no. of
%                               signals. The first two columns contain 2
%                               channels of measured EMG-signals and detected r peaks. Column 3
%                               contains the measured airway pressure, 4 the time
%                               and 5 the subject number.
%           fs              ->  Sampling rate of all signals
%           data            ->  Cell-matrix with same structure as
%                               'dataWithRPeaks' whose SNR improvement
%                               shall be calculated.
%   example:    analyzeSNRimprovement(dataWithRPeaks, 1024, ...
%                   dataKalmanFilter1, ...
%                   dataKalmanFilter2)
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

    % Number of cleaned data sets to analyze
    n = nargin-2;
    
    % skip first and last 10 seconds of analysis because these often contain artifacts
    skip_seconds = 10;    

    disp('Input...')
    SNRinput = calcEMGtoECGratio(dataWithRPeaks, fs, 0, [], skip_seconds);
    disp('First algo...')
    SNRoutput = calcEMGtoECGratio(dataWithRPeaks, fs, 0, data, skip_seconds);

    % Transform to dB
    SNRinput = 20 * log10(SNRinput);
    SNRoutput = 20 * log10(SNRoutput);
    
    % Throw away corrupted or otherwise unusable data sets
    SNRinput = SNRinput(:);
    SNRoutput = SNRoutput(:);

    % calculate actual SNR improvement
    SNRimprovement = zeros(length(SNRinput), n);
    SNRimprovement(:, 1) = SNRoutput - SNRinput;
    
    for ii = 2:n
        disp('Next algo...' )
        SNRoutput = calcEMGtoECGratio(dataWithRPeaks, fs, 0, varargin{ii - 1}, skip_seconds);
        SNRoutput = 20 * log10(SNRoutput);
        SNRoutput = SNRoutput(:);
        improvement = SNRoutput - SNRinput;
        SNRimprovement(:, ii) = improvement;
    end
end
