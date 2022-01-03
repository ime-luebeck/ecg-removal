% This program performs the Empirical Mode Decomposition according to the paper 
%       "On the HHT, its problems, and some solutions", Reference: Rato, R. T., Ortigueira, M. D., and Batista, A. G., 
%       Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.
%
%   Original Authors: Raul Rato (rtr@uninova.DOT.pt) and Manuel Ortigueira (mdortigueira@uninova.pt or mdo@fct.unl.pt) 
%
%   Modified by members of the IME, Universität zu Lübeck, Germany, January 2018:  Modified to only calculate the first IMF.
%
%--------------------------------------------------------------------------
%
%rParabEmd__L: Emd parabolic decomposition with extrapolated extrema
%v2.00
%
%   Usage:  rParabEmd= rParabEmd__L(x,qResol, qResid, qAlfa);
%           signal - input signal - must be a real vector
%           qResol - Resolution (in DBs: 10*log(WSignal/Bias energy))- normally between 40 and 60 dB 
%           qResid - Residual energy (in DBs: 10*log (WSignal/WqResidual))- normally between 40 and 60 dB
%           qAlfa  - Gradient step size (normally is set to 1)
%           rParabEmd    - relation matrix of IMF modes  (each as a line)
%                   with residual in last line.
%
%   Limitations:    NaN is not trapped
%
%   History:    V1.00 First version
%               V1.01 Count mismatch detection (Line 44) increased from 1 to 2
%               V2.00 time optimated          
%
%
% Original version: Copyright (c) 2008, Manuel Ortigueira
% All rights reserved.
% 
% WARNING: This software is a result of our research work and is supplied without any garanties.
%           We would like to receive comments on the results and report on bugs.
%
%           /* NoSPAM: Replace .DOT. with a dot (.) */
%                   (c) LaPAS-2007
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% - Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% - Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% - Neither the name of UNINOVA nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Modified version: Copyright 2019, Institute for Electrical Engineering in Medicine, 
% University of Luebeck
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

function signal_IMF = EMD(signal, qResol, qAlfa)

dbstop if warning
% error messages by using more or less than three input parameters
if(nargin~=3), error('EMG_single_IMF: Use with 3 inputs.'), end
if(nargout>1), error('EMG_single_IMF: Use with just one output.'), end
ArgCheck_s(signal, qResol, qAlfa)

% Actual computation -------------------------------------
signal_copy = signal(:);                  % get copy of the input signal
% loop to decompose the input signal into successive IMFs

% Originally, the following was contained in a while loop to extract the
% individual IMFs. We only want the first IMF, so no loop here.
signal_IMF = signal_copy; % at the beginning of the sifting process, kImf is the signal
IMF_length = length(signal_IMF);      % Signal length
IMF_Max= rGetPMaxs_s(signal_IMF);     % IMF_Max_Extrapol = [xM(M), yM(M)];
IMF_Min= rGetPMins_s(signal_IMF);     % IMF_Min_Extrapol = [xm(m), ym(m)];
IMF_Max_Extrapol= rPMaxExtrapol_s(IMF_Max, IMF_Min, IMF_length);
IMF_Min_Extrapol= rPMinExtrapol_s(IMF_Max, IMF_Min, IMF_length);
IMF_Max_Extrapol_length= length(IMF_Max_Extrapol);   
IMF_Min_Extrapol_length= length(IMF_Min_Extrapol);

% Check if IMF_Max_Extrapol_length and IMF_Min_Extrapol_length have the
% same length
if (abs(IMF_Max_Extrapol_length-IMF_Min_Extrapol_length)>2)
    disp('Debug: Max-Min count mismatch.');
end
if(sum(abs(diff(sign(IMF_Max_Extrapol(1:min(IMF_Max_Extrapol_length,IMF_Min_Extrapol_length),1)- IMF_Min_Extrapol(1:min(IMF_Max_Extrapol_length,IMF_Min_Extrapol_length),1)))))>0)
    disp('Debug: Max-Min sequence mismatch.');
end
if(sum(abs(diff(sign(IMF_Min_Extrapol(1:min(IMF_Max_Extrapol_length,IMF_Min_Extrapol_length),1)- IMF_Max_Extrapol(1:min(IMF_Max_Extrapol_length,IMF_Min_Extrapol_length),1)))))>0)
    disp('Debug: Max-Min reverse sequence mismatch.');
end
% interpolated vector IMF_M.._Interpol with values which are determined by cubic spline interpolation
IMF_Max_Interpol= spline(IMF_Max_Extrapol(:,1), IMF_Max_Extrapol(:,2), 1:IMF_length);          % Top envelop: IMF_Max_Interpol[n];
IMF_Min_Interpol= spline(IMF_Min_Extrapol(:,1), IMF_Min_Extrapol(:,2), 1:IMF_length);          % Down envelop: IMF_Min_Interpol[n];
IMF_Bias= (IMF_Max_Interpol+IMF_Min_Interpol)/2;               %  first bias estimate
while true(1)             % inner loop to find each IMF
    % calculate current IMF energy and bias energy
    IMF_energy = (signal_IMF)'*signal_IMF;                
    IMF_Bias_energy = IMF_Bias*IMF_Bias';                  
    if IMF_Bias_energy * IMF_energy < 0 
        warning('rParabEmd__L: Ooops, negative energy detected.');
    end
    if IMF_Bias_energy > 0 
        qResol_Temp = 10*log10(IMF_energy/IMF_Bias_energy); 
    else
        qResol_Temp = Inf;
    end
    if (qResol_Temp > qResol)
        %Resolution reached
        break; 
    end 
    %Resolution not reached. More work is needed
    signal_IMF = signal_IMF - qAlfa * IMF_Bias';                % subtract qAlfa bias from signal_IMF
    
    % work in next loop with signal_IMF minus qAlfa * IMF_Bias'
    IMF_Max = rGetPMaxs_s(signal_IMF);     % IMF_Max_Extrapol = [xM(M), yM(M)];
    IMF_Min = rGetPMins_s(signal_IMF);     % IMF_Min_Extrapol = [xm(m), ym(m)];
    IMF_Max_Extrapol = rPMaxExtrapol_s(IMF_Max, IMF_Min, IMF_length);
    IMF_Min_Extrapol = rPMinExtrapol_s(IMF_Max, IMF_Min, IMF_length);
    IMF_Max_Interpol = spline(IMF_Max_Extrapol(:,1), IMF_Max_Extrapol(:,2), 1:IMF_length);          % Top envelop: IMF_Max_Interpol[n];
    IMF_Min_Interpol = spline(IMF_Min_Extrapol(:,1), IMF_Min_Extrapol(:,2), 1:IMF_length);          % Down envelop: IMF_Min_Interpol[n];
    % new bias estimate
    IMF_Bias= (IMF_Max_Interpol+IMF_Min_Interpol)/2;              
end % Wend true
%
% IMF extraction + residual calculations omitted here.
% Outer while loop in original code stopped here.
% Some cleanup code for final IMF omitted here.
end 



function ArgCheck_s(signal, qResol, qAlfa)
	% Check all input parameters for correctness in terms of processing

	[signal_x, signal_y] = size(signal);
	% Check if input is an 1-dim vector
	if ((signal_x*signal_y) ~= max(signal_x,signal_y))
		error('rParabEmd__L: Input signal must be a one dim vector.');
	end
	if ((signal_x*signal_y) <= 1)
		error('rParabEmd__L: Input signal must be a vector.'); 
	end

	[signal_x,signal_y] = size(qResol);
	% Check if resolution is scalar and a strictly positive value
	if ( ~((signal_x == 1) & (signal_y == 1)) ) 
		error('rParabEmd__L: Input resolution must be a scalar.'); 
	end
	if ( qResol <= 0 ) 
		error('rParabEmd__L: Input resolution must strictly positive.'); 
	end

	[signal_x,signal_y] = size(qAlfa);
	% Check if Gradient step size qAlfa is a strictly positive scalar
	if ( ~((signal_x == 1) & (signal_y == 1)) ) 
		error('rParabEmd__L: qAlfa step must be a scalar.'); 
	end
	if ( qAlfa <= 0 ) 
		error('rParabEmd__L: qAlfa step  must be strictly positive.');  
	end
	end

	function rPMaxExtrapol= rPMaxExtrapol_s(IMF_Max, IMF_Min, IMF_length)
	% Time-mirrored top extrema (Parabolic Maxs) extrapolation

	%Init ------------------------------------
	IMF_Max = sortrows(IMF_Max); % assumes nothing on IMF_Max sort order
	IMF_Min = sortrows(IMF_Min); % assumes nothing on IMF_Min sort order

	IMF_Max_TopTime = IMF_Max(:,1); 
	IMF_Max_TopValue = IMF_Max(:,2);
	IMF_Min_TopTime = IMF_Min(:,1); 

	%Start extrapolation ---------------------
	if ( (IMF_Max_TopTime(1) == 1) && (IMF_Min_TopTime(1) == 1) )   
		disp ('rPMaxExtrapol_s: Poliextrema at signal''s start');
	elseif ( (IMF_Max_TopTime(1) < 1) || (IMF_Min_TopTime(1) < 1) )   
		disp ('rPMaxExtrapol_s: Invalid extrema at signal''s start');
	else
		IMF_Max_TopTime = [2 - IMF_Min_TopTime(1); IMF_Max_TopTime];     % New first Top at the (one based) specular Min
		IMF_Max_TopValue = [IMF_Max_TopValue(1); IMF_Max_TopValue];      % Same Val as old first Top
	end

	% End extrapolation -----------------------
	if ( (IMF_Max_TopTime(end) == IMF_length) && (IMF_Min_TopTime(end) == IMF_length) )   
		disp ('rPMaxExtrapol_s: Poliextrema at signal''s end');
	elseif ( (IMF_Max_TopTime(end ) > IMF_length) || (IMF_Min_TopTime(end) > IMF_length) )   
		disp ('rPMaxExtrapol_s: Invalid extrema at signal''s end');
	else
		IMF_Max_TopTime = [IMF_Max_TopTime; (2*IMF_length - IMF_Min_TopTime(end))];     % New last Top at the specular Min
		IMF_Max_TopValue = [ IMF_Max_TopValue; IMF_Max_TopValue(end)];                  % Same Val as old last Top 
	end

	rPMaxExtrapol = sortrows([IMF_Max_TopTime, IMF_Max_TopValue]); 

end



function rPMinExtrapol= rPMinExtrapol_s(IMF_Max, IMF_Min, IMF_length)
	% Time-mirrored down extrema (Parabolic Mins) extrapolation

	%Init ------------------------------------
	IMF_Max = sortrows(IMF_Max); %assumes nothing on rPM sort order
	IMF_Min = sortrows(IMF_Min); %assumes nothing on rPm sort order

	IMF_Max_TopTime = IMF_Max(:,1); 
	IMF_Min_TopTime = IMF_Min(:,1); 
	IMF_Min_TopValue = IMF_Min(:,2);

	%Start extrapolation ---------------------
	if ( (IMF_Max_TopTime(1) == 1) && (IMF_Min_TopTime(1) == 1) )
		disp ('rPMinExtrapol_s: Poliextrema at signal''s start');
	elseif ( (IMF_Max_TopTime(1) < 1) || (IMF_Min_TopTime(1) < 1) )
		disp ('rPMinExtrapol_s: Invalid extrema at signal''s start');
	else
		IMF_Min_TopTime = [2 - IMF_Max_TopTime(1); IMF_Min_TopTime];     % New first Dwn at the (one based) specular Max
		IMF_Min_TopValue = [IMF_Min_TopValue(1); IMF_Min_TopValue];      % Same Val as old first Dwn
	end

	% End extrapolation -----------------------
	if ( (IMF_Max_TopTime(end) == IMF_length) && (IMF_Min_TopTime(end) == IMF_length) )
		disp ('rPMinExtrapol_s: Poliextrema at signal''s end');
	elseif ( (IMF_Max_TopTime(end) > IMF_length) || (IMF_Min_TopTime(end) > IMF_length) )
		disp ('rPMinExtrapol_s: Invalid extrema at signal''s end');
	else
		IMF_Min_TopTime = [IMF_Min_TopTime; (2*IMF_length - IMF_Max_TopTime(end))];     % New last Dwn at the specular Max
		IMF_Min_TopValue = [ IMF_Min_TopValue; IMF_Min_TopValue(end)];          % Same Val as old last Dwn
	end

	rPMinExtrapol = sortrows([IMF_Min_TopTime, IMF_Min_TopValue]);

end



function rPMax = rGetPMaxs_s(signal)         
	%Get Parabolic Maxs, plateaus out

	signal_copy = signal(:);
	signal_length = length(signal_copy); 

	% Extract the peaks of input signal and determine their count
	[~, signal_peaks] = findpeaks(signal_copy);
	signal_peaks_count = length(signal_peaks);

	% Now we have the Maxs, lets get the Parabolic Maxs
	Xvalue_old = -Inf; 
	Yvalue_old = -Inf;
	signal_range = max(signal_copy) - min(signal_copy);

	% Initialise the arrays 
	[signal_time] = zeros(length(signal_copy),1);
	[signal_value] = zeros(length(signal_copy),1);
	Count = 0;
	% for all Maxs
	for k = 1 : signal_peaks_count    
		%sample_point_pre = -1; sample_point = 0; sample_point_after = 1;
		sample_point_pre = signal_copy(signal_peaks(k)-1);  % Sample point before
		sample_point = signal_copy(signal_peaks(k));    % Sample point, == signal_peaks(k)
		sample_point_after = signal_copy(signal_peaks(k)+1);  % Sample point after
		D = (-4 * sample_point + 2*sample_point_pre + 2*sample_point_after);
		if (D == 0) 
			Xvalue = signal_peaks(k);
		else
			Xvalue = signal_peaks(k) + (sample_point_pre - sample_point_after) / D; 
		end 
		D = (-16*sample_point+ 8*sample_point_pre+ 8*sample_point_after);
		if (D == 0) 
			Yvalue = sample_point;
		else
			Yvalue = sample_point + (2*sample_point_after*sample_point_pre - sample_point_pre*sample_point_pre - sample_point_after*sample_point_after) / D; 
		end
		% Lets check for double maxima
		if ( (Xvalue == Xvalue_old)||(abs(Yvalue - Yvalue_old)/abs(Xvalue - Xvalue_old)) > (2*signal_range))       
			Xvalue = (Xvalue + Xvalue_old) / 2; 
			Yvalue = max(Yvalue,Yvalue_old);   %Double found
			signal_time(Count) = Xvalue;
			signal_value(Count) = Yvalue;
			Count = Count + 1;
		else
			Count = Count + 1;
			signal_time(Count,1) = Xvalue;
			signal_value(Count,1) = Yvalue;
		end 
		Xvalue_old = Xvalue;
		Yvalue_old = Yvalue;
	end % for k = 1:signal_peaks_count

	[signal_time] = signal_time(1 : Count-1,1);
	[signal_value] = signal_value(1 : Count-1,1);

	if signal_peaks_count > 0
		if ( signal_copy(1) >= signal_value(1) )
			signal_time = [1; signal_time];  
			signal_value =[signal_copy(1); signal_value ];    %Start must be included as a Max
		end
		if ( signal_copy(end) >= signal_value(end))
			signal_time = [signal_time; signal_length];  
			signal_value = [signal_value; signal_copy(end)];   %End must be included as a Max
		end
	end

	if signal_peaks_count == 0
		if ( signal_copy(1) > signal_copy(2) )
			signal_time = [1; signal_time];  
			signal_value = [signal_copy(1); signal_value];    %Start must be included as a Max
		end
		if ( signal_copy(end) > signal_copy(end-1))
			signal_time = [signal_time; signal_length];  
			signal_value = [signal_value; signal_copy(end)];   %End must be included as a Max
		end
	end
	if signal_peaks_count < 0
		error('rGetPMaxs_s: Invalid MaxCnt value');
	end


	rPMax = sortrows([signal_time, signal_value]);
end


 
function rPMin = rGetPMins_s(signal)         
	%Get Parabolic Mins, plateaus out

	signal_copy = signal(:);
	signal_length = length(signal_copy); 

	[~, signal_peaks] = findpeaks(-signal_copy);
	signal_peaks_count = length(signal_peaks);

	% Now we have the Mins, lets get the Parabolic Mins
	Xvalue_old = -Inf; 
	Yvalue_old = -Inf;
	signal_range = max(signal_copy) - min(signal_copy);

	% Initialise arrays
	[signal_time] = zeros(length(signal_copy),1);
	[signal_value] = zeros(length(signal_copy),1);
	Zaehlvar = 0;

	for k = 1 : signal_peaks_count     %for all Mins
		%sample_point_pre = -1; sample_point = 0; sample_point_after = 1;
		sample_point_pre = signal_copy(signal_peaks(k)-1);  % Sample point before
		sample_point = signal_copy(signal_peaks(k));    % Sample point, == signal_peaks(k)
		sample_point_after = signal_copy(signal_peaks(k)+1);  % Sample point after
		D = (-4 * sample_point + 2*sample_point_pre + 2*sample_point_after);
		if (D == 0)
			Xvalue = signal_peaks(k);
		else
			Xvalue = signal_peaks(k) + (sample_point_pre-sample_point_after) / D; 
		end 
		D = (-16 * sample_point + 8*sample_point_pre + 8*sample_point_after);
		if (D == 0) 
			Yvalue = sample_point;
		else
			Yvalue = sample_point + (2*sample_point_after*sample_point_pre - sample_point_pre*sample_point_pre - sample_point_after*sample_point_after)/D; 
		end
		% Lets check for double minima
		if ( (Xvalue == Xvalue_old)||(abs(Yvalue - Yvalue_old)/abs(Xvalue - Xvalue_old)) > (2*signal_range) )     
			Xvalue = (Xvalue + Xvalue_old)/2; 
			Yvalue = min(Yvalue,Yvalue_old);   % Double found
			signal_time(Zaehlvar) = Xvalue;
			signal_value(Zaehlvar) = Yvalue;
			Zaehlvar = Zaehlvar + 1;
		else
			Zaehlvar = Zaehlvar + 1;
			signal_time(Zaehlvar,1) = Xvalue;
			signal_value(Zaehlvar,1) = Yvalue;

		end 
		Xvalue_old = Xvalue; 
		Yvalue_old = Yvalue;
	end % for k = 1:signal_peaks_count

	[signal_time] = signal_time(1:Zaehlvar-1,1);
	[signal_value] = signal_value(1:Zaehlvar-1,1);

	if signal_peaks_count > 0
		if ( signal_copy(1) <= signal_value(1) )
			signal_time = [1; signal_time];  
			signal_value = [signal_copy(1); signal_value ];    %Start must be included as a Min
		end
		if ( signal_copy(end) <= signal_value(end))
			signal_time = [signal_time; signal_length];  
			signal_value = [signal_value; signal_copy(end)];   %End must be included as a Min
		end
	end

	if signal_peaks_count == 0
		if ( signal_copy(1) < signal_copy(2) )
			signal_time = [1; signal_time];  
			signal_value = [signal_copy(1); signal_value];    %Start must be included as a Min
		end
		if ( signal_copy(end) < signal_copy(end-1))
			signal_time = [signal_time; signal_length];  
			signal_value = [signal_value; signal_copy(end)];   %End must be included as a Min
		end
	end
	if signal_peaks_count < 0
		error('rGetPMins_s: Invalid MinCnt value');
	end


	rPMin = sortrows([signal_time, signal_value]);
end