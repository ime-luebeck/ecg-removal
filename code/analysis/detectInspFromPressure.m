function activity = detectInspFromPressure(pressure, fs)
%DETECTINSPFROMPRESSURE Detect inspirations in pressure signal
%
%   Inspiration (negative pressure) defines EMG activity. A safety phase of
%   0.05 seconds is added. Output of this function is a mask array.
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

    safetyPhase = ceil(0.05 * fs);

    thresh = -0.45;
    zeropm = find(diff(sign(pressure - thresh)) < 0);
    zeromp = find(diff(sign(pressure - thresh)) > 0);
    if zeromp(1) < zeropm(1)
        zeromp = zeromp(2:end);
    end
    if zeropm(end) > zeromp(end)
        zeropm = zeropm(1:end-1);
    end

    M = length(zeropm);
    activity = zeros(size(pressure));
    for i = 1:M
        activity(zeropm(i):zeromp(i) + safetyPhase) = 1;
    end

end