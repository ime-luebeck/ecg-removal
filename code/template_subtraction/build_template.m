function [template] = build_template(rWaveIndex, lengthLeft, lengthRight, ...
    signal, num , i, centered)
%BUILD_TEMPLATE Basic template formation by simple superposition
%   Returns an ECG artifact template that is formed by averaging
%   num ecg occurences, centered around rWaveIndex and of length lengthLeft
%   + lengthRight + 1. No scaling, segmentation, etc. done.
%   ECG beats to use are chosen next to beat i.
%   By default the template is calculated in a centered fashion, i.e., half
%   of the num beats are to the right and half to the left of beat i, and 
%   beat i is also included in the template formation.
%   To obtain causal template formation, pass centered = 'false'.
%   In that case, num beats to the left of beat i are used.
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

if nargin < 7 || isempty(centered)
    centered = true;
end

if centered
    beat_idx_start = i - ceil(num/2);
    beat_idx_stop = beat_idx_start + num - 1;
else
    beat_idx_start = i - num;
    beat_idx_stop = i - 1;
end

[m, n] = size(signal);
if m > n
    signal = signal';
    was_col_vec = true;
else
    was_col_vec = false;
end

template = zeros(1, lengthLeft + lengthRight + 1);
beatcount = zeros(size(template));

for beat_idx = beat_idx_start:beat_idx_stop
    if beat_idx < 1 || beat_idx > length(rWaveIndex)
        continue; % no data available
    else
        leftIdx = rWaveIndex(beat_idx) - lengthLeft;
        rightIdx = rWaveIndex(beat_idx) + lengthRight;
        if leftIdx < 1
            mask = 1-leftIdx+1:length(template);
            template(mask) = template(mask) + signal(1:rightIdx);
        elseif rightIdx <= length(signal)
            mask = 1:length(template);
            template = template + signal(leftIdx:rightIdx);
        else
            mask = 1:length(template)-(rightIdx-length(signal));
            template(mask) = template(mask) + signal(leftIdx:length(signal));
        end
    end
    
    % Count the number of beats that has contributed to the individual
    % samples in the template. This number may be different for different
    % samples due to incomplete beats at signal start and end.
    beatcount(mask) = beatcount(mask) + 1;
end

% Calculate average template
template = template ./ beatcount;
template(beatcount==0) = 0;

assert(all(size(template) == [1, lengthLeft+lengthRight+1]));
if was_col_vec
    template = template';
end

end

