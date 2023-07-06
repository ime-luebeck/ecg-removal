function varargout = filter_signals(sigIn, algo, channels, varargin)
%FILTERSIGNALS Run list of filters on set of signals
%   Runs a filter algorithm on a cell array of signals.
%
%   INPUT ARGUMENTS
%
%   sigIn: Cell array (nxm) of n datasets, each consisting of m measurement
%          channels (e.g. EMG1, EMG2, Paw).
%
%   algo: Function handle to the filtering algorithm to apply.
%
%   channels: List of channels (i.e., columns of sigIn) to apply the filtering 
%             algorithm to.
%
%   passAuxPerChannel: (Optional, positional) If this is nonzero, the signal
%                      sigIn{ii, chan}(:, passAuxPerChannel) is passed as an 
%                      additional input vector to the filter function. If 0 
%                      (default), no additional vector is passed.
%
%   passAuxPerRecording: (Optional, positional) If this is nonzero, the
%                        signal sigIn{ii, passAuxPerRecording}):, 1) is passed
%                        as an additional input vector to the filter function.
%                        If 0 (default), no additional vector is passed.
%
%   OUTPUTS
%
%   sigOut: Cell array of size (n x length(channels)) that contains the
%           filtered signals.
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

if nargin > 3
    passAuxPerChannel = varargin{1};
else
    passAuxPerChannel = 0;
end

if nargin > 4
    passAuxPerRecording = varargin{2};
    if passAuxPerRecording ~=0
        assert(passAuxPerChannel == 0);
    end
else
    passAuxPerRecording = 0;
end

n = size(sigIn, 1);


% Some algorithms return multiple results that should all be passed on.
% "nargout" yields the number of output arguments this function was called
% with. We expect to collect the same amount of output arguments from each
% call to "algo".
sigsOut = cell(nargout, 1);
for ii = 1:nargout
    sigsOut{ii} = sigIn;
end

% Loop over data sets
for ii = 1:n
    % Loop over channels in each data set that the filtering shall be applied to
    disp(['Dataset ', num2str(ii)]);
    for jj = 1:length(channels)
        chan = channels(jj);
        sigsOutLoc = cell(nargout, 1);
        
		if passAuxPerChannel == 0
			if passAuxPerRecording == 0
				[sigsOutLoc{:}] = algo(sigIn{ii, chan}(:, 1));
            else
				[sigsOutLoc{:}] = algo(sigIn{ii, chan}(:, 1), sigIn{ii, passAuxPerRecording}(:, 1));
			end
		else
			[sigsOutLoc{:}] = algo(sigIn{ii, chan}(:, 1), sigIn{ii, chan}(:, passAuxPerChannel));
		end

        for kk = 1:nargout
            sigsOut{kk}{ii, chan} = sigsOutLoc{kk};
        end
            
    end
end

varargout = sigsOut;