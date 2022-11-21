function ECGremovedSignal = template_subtraction(signal, rPeaks, varargin)
%TEMPLATE_SUBTRACTION Basic implementation with L2 error minimization.
%   Positions whole template by correlation with measurements. Template
%   amplitude is adjusted by a linear scaling factor to minimize L2 error
%   between adjusted template and measurements. This is done by means of a
%   simple linear regression between signal and template, which also
%   accounts for constant offset in either signal. Template width is
%   restricted to width of current RR-interval by simply cutting off the
%   remaining samples at the edges.
%
% The only optional argument is the number of observed QRS complexes to 
% build the current template from, default: 40.
%
% Notice that this implementation is completely agnostic to the sampling
% rate.
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

if nargin > 2
    num = varargin{1};
else
    num = 40;
end

rPeakIdces = find(rPeaks);

signal = signal(:);
N = size(signal);
numQRS = length(rPeakIdces);

subtractionTemplate = zeros(N);

% distance between R waves
distance = rPeakIdces(2:end) - rPeakIdces(1:end-1);

% center between QRS complexes (row vector)
center = floor(distance/2) + rPeakIdces(1:end-1);
center = [1; center(:); length(signal)];

% averaged template for every QRS complex
for i = 1:numQRS
    lengthLeft = rPeakIdces(i)-center(i);
    lengthRight = center(i+1)-rPeakIdces(i);
    curr_beat = signal(rPeakIdces(i)-lengthLeft:rPeakIdces(i)+lengthRight);
    % Build current template, using signal information including current artefact
    currentTemplate = build_template(rPeakIdces, lengthLeft, lengthRight, signal, num, i);

    % Position template to maximize correlation between template and current signal segment.
    [correlation, lags] = xcorr(curr_beat,currentTemplate, 'coeff');
    [~, indexmax] = max(correlation);
    lag = lags(indexmax); % lag > 0 => template must be shifted to the right
    
    if lag > lengthLeft || lag < -lengthRight
        lag = 0;
    end
    
    if lag > 0  % Cut off right end and move it to the left side
        rightEnd = currentTemplate(end-lag+1:end);
        currentTemplate(lag+1:end) = currentTemplate(1:end-lag);
        currentTemplate(1:lag) = rightEnd;
    elseif lag < 0 % Cut off left end and move it to the right side
        leftEnd = currentTemplate(1:-lag);
        currentTemplate(1:end+lag) = currentTemplate(-lag+1:end);
        currentTemplate(end+lag+1:end) = leftEnd;
    end
   
    % Scale the current template by means of linear regression: (signal - a * template - b).^2 != min
    coeffs = [ones(length(currentTemplate), 1) currentTemplate] \ curr_beat;
    currentTemplate = coeffs(2) * currentTemplate;    

    subtractionTemplate(rPeakIdces(i)-lengthLeft:rPeakIdces(i)+lengthRight) = currentTemplate;
end

ECGremovedSignal = signal - subtractionTemplate;
end

