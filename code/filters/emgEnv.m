function sigOut = emgEnv(sigIn, length, type)
%EMGENV Calculates Envelope of EMG signal by specified method.
%   INPUT:  sigIn   ->  Vector with EMG signal.
%           length  ->  windowlength for e.g. Moving Mean
%           type    ->  Possible methods are: 
%                       'rms', root-mean-squares
%                       'med', Median
%                       'abs', Moving Mean of absolute signal
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

    if type == 'rms'
        sigOut = sqrt(movmean(sigIn.^2, length));
    elseif type == 'med'
        sigOut = sqrt(movmedian(sigIn.^2, length));
    elseif type == 'abs'
        sigOut = movmean(abs(sigIn), length);
    elseif type == 'cwt'
        sigOut = cwtEmgEnvelope(sigIn, length);
    end
end