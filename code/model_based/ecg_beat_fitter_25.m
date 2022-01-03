function [OptimumParams, Model, P_std, C_std, T_std] = ecg_beat_fitter_25(meanECG, meanPhase, varargin)
%ECG_BEAT_FITTER_25 Parameter estimation of synthetic ECG Model with
%separation of P-wave, QRS-complex and T-wave sections.
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

if nargin == 2
    [Pmean, Cmean, Tmean] = ecg_windowing(meanECG, meanPhase);
elseif nargin == 5
    a = varargin{1};
    b = varargin{2};
    gamma = varargin{3};
    [Pmean, Cmean, Tmean] = ecg_windowing(meanECG, meanPhase, a, b, gamma);
end

[~,PInd] = max(Pmean);

Cmean2 = movmean(Cmean, 10, 'omitnan');
CInd = sort_peaks(Cmean2);
n = min(length(CInd), 3);
CInd = CInd(1:n);

[~,TInd] = max(Tmean);

Ind = [max(PInd-25, 1), PInd, CInd(:)', TInd, min(TInd+25, length(meanECG))];
tetai = meanPhase(Ind);
alphai = meanECG(Ind);
bi = .04 * ones(size(alphai));

m = length(Ind);
OptimumParams = zeros(3*m, 1);

options = optimset('TolX', 1e-4, 'TolFun', 1e-4, 'MaxIter', 100);

P_initParams = [alphai(1:2) bi(1:2) tetai(1:2)];
P_optimumParams = nlinfit(meanPhase', Pmean, @ecg_model, P_initParams, options);
OptimumParams([1:2, m+1:m+2, 2*m+1:2*m+2]) = P_optimumParams;
P_std = ecg_model_error(P_optimumParams, Pmean, meanPhase, 0);

C_initParams = [alphai(3:end-2) bi(3:end-2) tetai(3:end-2)];
C_optimumParams = nlinfit(meanPhase', Cmean, @ecg_model, C_initParams, options);
OptimumParams([3:m-2, m+3:2*m-2, 2*m+3:3*m-2]) = C_optimumParams;
C_std = ecg_model_error(C_optimumParams, Cmean, meanPhase, 0);

T_initParams = [alphai(end-1:end) bi(end-1:end) tetai(end-1:end)];
T_optimumParams = nlinfit(meanPhase', Tmean, @ecg_model, T_initParams, options);
OptimumParams([m-1:m, 2*m-1:2*m, 3*m-1:3*m]) = T_optimumParams;
T_std = ecg_model_error(T_optimumParams, Tmean, meanPhase, 0);

Model = ecg_model_error(OptimumParams, meanECG, meanPhase, 1);

end