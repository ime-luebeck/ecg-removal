function plot_results(dataWithRPeaks, env_method, filename, varargin)
%PLOT_RESULTS Plots multiple measured signals (structure of
%'dataWithRPeaks').
%   INPUT:  dataWithRPeaks  ->  nx5 cell-matrix, where n is the no. of
%                               signals. The first two columns contain 2
%                               channels of measured EMG-signals and
%                               detected r peaks. Column 3 contains the
%                               measured airway pressure (can be empty), 
%                               4 the time (ms) and 5 the subject number.
%           env_method      ->  integer 0 to 3. If 0 the no envelope is
%                               plotted. 1 to 3 are different methods of
%                               calculating the envelope of the ECG-removed
%                               signals.
%           filename        ->  filename to export plots to
%           algo name       ->  String with description of 'data'.
%           data            ->  Cell-matrix with same structure as
%                               'dataWithRPeaks' with signals that shall
%                               be plotted as well.
%   example:    plotResults(dataWithRPeaks, 0, ...
%                   'TS', dataTS, ...
%                   'EKS2', dataEKS2)
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

n = 100;
n_data = floor((nargin - 3) / 2);
mode = [];
n_rows = n_data + 1;
n_cols = 2;

if n_data < 1
    mode = 'measurement';
    env_method = 0;
    n_rows = 3;
    n_cols = 1;
end

slider_size = 20;
CONST.length = 0.8 / n_cols;
CONST.lmargin = (1 - n_cols * CONST.length) / (n_cols + 1);

for subject_id = 1:size(dataWithRPeaks, 1)   % iterate over subjects
    
    fig = figure('name', ['subject: ', dataWithRPeaks{subject_id,5}]);
    time = dataWithRPeaks{subject_id, 4}(n:end) / 1000; % time [seconds]
    
    axes_handles = cell(n_rows, n_cols);
    for ax = 1:numel(axes_handles)
        axes_handles{ax} = axes;
    end
    
    if strcmp(mode, 'measurement')        
        plot(axes_handles{1}, time, dataWithRPeaks{subject_id, 3}(n:end, 1))
        ylabel(axes_handles{1}, 'Pressure (kPa)')
        
        plot(axes_handles{2}, time, dataWithRPeaks{subject_id, 1}(n:end, 1))
        ylabel(axes_handles{2}, 'EMG 1 (mV)')
        
        plot(axes_handles{3}, time, dataWithRPeaks{subject_id, 2}(n:end, 1))
        ylabel(axes_handles{3}, 'EMG 2 (mV)')
    else
        % Plot raw EMG signals with base
        plot(axes_handles{1, 1}, time, dataWithRPeaks{subject_id, 1}(n:end, 1))
        ylabel(axes_handles{1, 1}, 'Raw EMG (mV)')
        title(axes_handles{1, 1}, 'Channel 1')
        
        plot(axes_handles{1, 2}, time, dataWithRPeaks{subject_id, 2}(n:end, 1))
        title(axes_handles{1, 2}, 'Channel 2')
        
        if ~isempty(env_method)
        % Calculate and plot envelope of the ECG-removed signals instead.
            for row = 2:n_rows
                for col = 1:n_cols
                    env = emgEnv(varargin{2 * (row - 1)}{subject_id, col}(n:end, 1), 100, env_method);
                    plot(axes_handles{row, col}, time, env, ...
                        'Color', [115, 115, 115] / 256);
                end
                algo_label = strcat('EMG envelope (', varargin{2 * (row - 1) -1}, ') (mV)');
                ylabel(axes_handles{row, 1}, algo_label)
            end
        else
            % Plot raw EMG signals after ECG removal
            for ii = 1:n_data
                jj = 2 * ii - 1;

                plot(axes_handles{ii + 1, 1}, time, ...
                    varargin{jj + 1}{subject_id, 1}(n:end, 1)) 

                plot(axes_handles{ii + 1, 2}, time, ...
                    varargin{jj + 1}{subject_id, 2}(n:end, 1))

                algo_label = strcat('EMG (', varargin{jj}, ') (mV)');
                ylabel(axes_handles{ii + 1, 1}, algo_label)
            end
        end
    end

    CONST.hmargin = 0.02;
    CONST.height = 1/n_rows - 2 * CONST.hmargin;
    for row = 1:n_rows
        for col = 1:n_cols
           axes_handles{row, col}.Position = [ ...
               col * CONST.lmargin + (col - 1) * CONST.length, ...
               CONST.hmargin + (n_rows - row) * ...
               (CONST.height + 2 * CONST.hmargin), ...
               CONST.length, CONST.height];            
           xlabel(axes_handles{row, col}, 'time (s)')
        end
    end
    linkaxes(findall(fig, 'type', 'axes'), 'x')
    
    % export figure
    fig_width = 10;
    fig_height = 2 * n_rows;
    print_fig_to_png(fig, [filename, '_subject', dataWithRPeaks{subject_id, 5}], fig_width, fig_height);
    xlim([500, 545]);
    % export figure zoomed in
    print_fig_to_png(fig, [filename, '_subject', dataWithRPeaks{subject_id, 5}, '_zoom'], fig_width, fig_height);
    xlim auto;
    
    % now (after export) switch to scrollable figure design, making some
    % rows invisible initially
    on_screen = min(4, n_rows);
    CONST.hmargin = 0.04;
    CONST.height = 1/on_screen - 2 * CONST.hmargin;
    for row = 1:n_rows
        for col = 1:n_cols
           axes_handles{row, col}.Position = [ ...
               col * CONST.lmargin + (col - 1) * CONST.length, ...
               CONST.hmargin + (on_screen - row) * ...
               (CONST.height + 2 * CONST.hmargin), ...
               CONST.length, CONST.height];
        end
    end
           
    % scrollable figure, if there are many subplots
    if on_screen < n_rows
        FigurePosition = fig.Position;
        SliderPositionX = FigurePosition(3) - slider_size;
        Slider = uicontrol(...
            'Style','slider',...
            'Position',[SliderPositionX, 0, slider_size, FigurePosition(4)],...
            'Min',- (n_rows - on_screen) * (CONST.height + 2*CONST.hmargin),...
            'Max', 0,...
            'Value', 0,...
            'SliderStep', [1 - on_screen/n_rows, 2 * on_screen/n_rows], ...
            'Callback',{@slidermove, axes_handles, on_screen, CONST});
        fig.SizeChangedFcn = {@refreshSlider, Slider, slider_size};
    end
end

% Figure Size Change Callback function
function refreshSlider(source, ~, Slider, SLIDER_SIZE)
    FigurePosition_ = source.Position;
    SliderPositionX_ = FigurePosition_(3) - SLIDER_SIZE;
    Slider.Position = [SliderPositionX_, 0, SLIDER_SIZE, FigurePosition_(4)];
end

% Slider Move Callback function
function slidermove(source, ~, axes_handles, on_screen, CONST)
    for rr = 1:size(axes_handles, 1)
       for cc = 1:size(axes_handles, 2)
          axes_handles{rr, cc}.Position(2) = CONST.hmargin + ...
              (on_screen - rr) * (2*CONST.hmargin + CONST.height) - ...
              source.Value; 
       end
    end
end

end

