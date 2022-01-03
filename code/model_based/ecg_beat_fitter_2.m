function [ OptParams, Model,num ] = ecg_beat_fitter_2( MeanECG, MeanPhase, num )
%ECG_BEAT_FITTER_2 Parameter estimation of synthetic ECG Model 
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

meanMeanECG = movmean(MeanECG, 10, 'omitnan');
[peaks,peakIndex] = findpeaks(abs(meanMeanECG));
peakmatrix = [peaks(:) peakIndex(:)];
sortedpeakmatrix = flipud(sortrows(peakmatrix));
sortedIndex = sortedpeakmatrix(:, 2);
num = min(length(sortedIndex), num);
I = sortedIndex(1:num);

tetai = MeanPhase(I);
alphai = 1.2 * MeanECG(I);
bi = .04 * ones(size(alphai));

options = optimset('TolX', 1e-4, 'TolFun', 1e-4, 'MaxIter', 100);
InitParams = [alphai bi tetai];
OptParams = nlinfit(MeanPhase, MeanECG, @ecg_model, InitParams, options);
Model = ecg_model_error(OptParams, MeanECG, MeanPhase, 1);

end

