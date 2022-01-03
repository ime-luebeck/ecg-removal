function print_fig_to_png(fig, filename, x_width_in, y_width_in)

    set(fig, 'PaperUnits', 'inches');
    if nargin < 4
        % Golden Ratio for plotting
        GR=1.618;
        y_width_in = x_width_in / GR;
    end
    set(fig, 'PaperPosition', [0 0 x_width_in y_width_in]);
    print(fig, filename, '-dpng');