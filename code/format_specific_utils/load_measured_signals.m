function data = load_measured_signals(directory)
%LOAD_MEASURED_SIGNALS Stores the selected data from the directory into the
%cell array 'data'.
%   INPUT
%       directory:  directory of the measured signals
%   OUTPUT
%       data =      {channel1, channel2, pressure, time,'01'     -> of testperson 01
%                    channel1, channel2, ...     ,     ,'02'}    -> of testperson ...
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

% load all measured signals
files = dir(directory);
fileIndex = find(~[files.isdir]);
n = length(fileIndex);
data = cell(n,5);
for i = 1:n
    fileName = files(fileIndex(i)).name;
    no = fileName(end-4);
    t = readtable(fullfile(directory, fileName));
    data{i,1} = t.Ch1_mV;
    data{i,2} = t.Ch2_mV;
    data{i,3} = t.Pressure_kPa;
    data{i,4} = t.Time_ms;
    data{i,5} = no;
end

