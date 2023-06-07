function ECGremovedSignal = adaptive_template_subtraction(signal, rPeaks, fs, varargin)
%ADAPTIVE_TEMPLATE_SUBTRACTION Advanced implementation that handles and
%   fits QRS complex and the remaining template separately. QRS section of
%   template (R peak +/- 40ms) is positioned by linear correlation, then
%   shrunk / stretched by pchip interpolation to fit the signal optimally.
%   To achieve this, QRS section width is increased/decreased successively
%   by one sample on each side, up to max. +/- 10 samples on each side
%   (i.e., 21 settings total). The resulting QRS section template is
%   calculated by (pchip) interpolation of the original QRS section
%   template. The (unmodified) remainder of the template on either side is
%   cut off to match the desired total signal length. Total squared error
%   is then calculated for each of these 21 differently stretched template
%   versions by correlation and simple linear regression as in the simpler 
%   implementations; the optimal solution is selected.
%
% The only optional argument is the number of observed QRS complexes to 
% build the current template from, default: 40.
%
% Copyright 2019 Institute for Electrical Engineering in Medicine, 
% University of Luebeck
% Eike Petersen
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


if nargin > 3
    num = varargin{1};
else
    num = 40;
end

rPeakIdces = find(rPeaks);

signal_sizes = size(signal);
if signal_sizes(1) == 1
    was_row = true;
else
    was_row = false;
end

signal = signal(:);
N = size(signal);
numQRS = length(rPeakIdces);

qrsIncMax = 10;

% Assumed QRS length is 110ms with R peak at center.
qrsLengthLeft = ceil(0.055 * fs);
qrsLengthRight = qrsLengthLeft;

subtractionTemplate = zeros(N);

% distance between R waves
distance = rPeakIdces(2:end) - rPeakIdces(1:end-1);

% center between QRS complexes (row vector)
center = floor(distance/2) + rPeakIdces(1:end-1);
center = [1; center(:); length(signal)];

optIdcesCounts = zeros(2*qrsIncMax+1, 1);

% averaged template for every QRS complex
for i = 1:numQRS
    lengthLeft = rPeakIdces(i)-center(i);
    lengthRight = center(i+1)-rPeakIdces(i);
    next = signal(rPeakIdces(i)-lengthLeft:rPeakIdces(i)+lengthRight);
    % Build current template, using signal information including current artefact
    currentTemplate = build_template(rPeakIdces, lengthLeft, lengthRight, signal, num, i);
    assert(lengthLeft + lengthRight + 1 == length(currentTemplate));
    
    qrsLeftEndIdx = lengthLeft + 1 - qrsLengthLeft;
    qrsRightEndIdx = lengthLeft + 1 + qrsLengthRight;
    if qrsLeftEndIdx > qrsIncMax && qrsRightEndIdx < length(currentTemplate) - qrsIncMax
        qrs = currentTemplate(qrsLeftEndIdx:qrsRightEndIdx);

        templates = zeros(length(currentTemplate), 2*qrsIncMax+1);
        L2error = zeros(2*qrsIncMax+1, 1);

        for qrsInc = -qrsIncMax:qrsIncMax
            % Construct modified template: Stretch QRS complex to length(qrs) + 2*qrsmod samples (by adding qrsmod 
            % samples on each side of the R peak)
            inc = (length(qrs) - 1) / (length(qrs) - 1 + 2*qrsInc);
            qrsMod = interp1(1:length(qrs), qrs, 1:inc:length(qrs), 'pchip');
            qrsMod = qrsMod';
            assert(length(qrsMod) == length(qrs) + 2 * qrsInc);

            % Concatenate with the remainder of the current template; cut off or add zeroes to reach target length
            if qrsInc > 0
                % QRS width has been increased; therefore we must cut off
                leftPart = currentTemplate(1+qrsInc:qrsLeftEndIdx-1);
                rightPart = currentTemplate(qrsRightEndIdx+1:end-qrsInc);
            else
                % QRS width has been reduced; therefore we must augment by zeroes
                leftPart = [zeros(-qrsInc, 1); currentTemplate(1:qrsLeftEndIdx-1)];
                rightPart = [currentTemplate(qrsRightEndIdx+1:end); zeros(-qrsInc, 1)];
            end
            currentTemplateMod = [leftPart; qrsMod; rightPart];
            assert(length(currentTemplateMod) == length(currentTemplate));
            
            % Position template to maximize correlation between template and current signal segment.
            [correlation, lags] = xcorr(next, currentTemplateMod, 'coeff');
            [~, indexmax] = max(correlation);
            lag = lags(indexmax); % lag > 0 => template must be shifted to the right

            if lag > lengthLeft || lag < -lengthRight
                lag = 0;
            end

            if lag > 0  % Cut off right end and move it to the left side
                rightEnd = currentTemplateMod(end-lag+1:end);
                currentTemplateMod(lag+1:end) = currentTemplateMod(1:end-lag);
                currentTemplateMod(1:lag) = rightEnd;
            elseif lag < 0 % Cut off left end and move it to the right side
                leftEnd = currentTemplateMod(1:-lag);
                currentTemplateMod(1:end+lag) = currentTemplateMod(-lag+1:end);
                currentTemplateMod(end+lag+1:end) = leftEnd;
            end
			lengthLeftMod = length(leftPart) + lag;
			lengthRightMod = length(rightPart) - lag;
			
            if lengthLeftMod >= 0 && lengthRightMod >= 0
                baseLine = zeros(length(currentTemplateMod), 1);
                % Scale the different signal sections separately by means of linear regression: 
                % sum((signal - a * template - b).^2) != min
                coeffsLeft = [ones(lengthLeftMod, 1) currentTemplateMod(1:lengthLeftMod)] \ next(1:lengthLeftMod);
                leftPartMod = coeffsLeft(2) * currentTemplateMod(1:lengthLeftMod);
                baseLine(1:lengthLeftMod) = coeffsLeft(1);

                qrsMask = lengthLeftMod+1:length(currentTemplateMod)-lengthRightMod;
                coeffsQrs = [ones(length(qrsMask), 1) currentTemplateMod(qrsMask)] \ next(qrsMask);
                qrsMod = coeffsQrs(2) * currentTemplateMod(qrsMask);
                baseLine(qrsMask) = coeffsQrs(1);

                coeffsRight = [ones(lengthRightMod, 1) currentTemplateMod(end-lengthRightMod+1:end)] \ next(end-lengthRightMod+1:end);
                rightPartMod = coeffsRight(2) * currentTemplateMod(end-lengthRightMod+1:end);
                baseLine(end-lengthRightMod+1:end) = coeffsRight(1);

                currentTemplateMod = [leftPartMod; qrsMod; rightPartMod];
            else
                % One of the signal segments is degenerated; just scale the whole signal at once.
                coeffs = [ones(length(currentTemplateMod), 1) currentTemplateMod] \ next;
                currentTemplateMod = coeffs(2) * currentTemplateMod;
                baseLine = ones(length(currentTemplateMod), 1) * coeffs(1);
            end
            
            % Save these results for later consideration
            templates(:, qrsInc + qrsIncMax + 1) = currentTemplateMod;
            L2error(qrsInc + qrsIncMax + 1) = sum((next - currentTemplateMod - baseLine).^2);
        end

        [~, optIdx] = min(L2error);
        optIdcesCounts(optIdx) = optIdcesCounts(optIdx) + 1;
        subtractionTemplate(rPeakIdces(i)-lengthLeft:rPeakIdces(i)+lengthRight) = templates(:, optIdx);
    else
        % Weirdly short segment; we're probably totally lost anyway. Just carry on...
        subtractionTemplate(rPeakIdces(i)-lengthLeft:rPeakIdces(i)+lengthRight) = currentTemplate;
    end
end

ECGremovedSignal = signal - subtractionTemplate;

if was_row
    ECGremovedSignal = ECGremovedSignal';
end

end

