#! /usr/bin/env python
"""
Provide a directory containing ModEm logfiles to have them
joined into a single file. This joint logfile is saved as
'{directory_name}.log'.

Logfiles are sorted alphanumerically before they are
joined, so ensure the files in the directory are named
correctly to achieve the desired ordering.

By default, the 'rms' value is plotted, but any value that
appears in the logfile per iteration as '{metric}={value}'
can be plotted'. Plot is saved as '{metrc}.png'.

See 'mtpy.modeling.modem.plot_rms_iterations' for implementation
and a command line interface.
"""
from mtpy.modeling.modem.plot_rms_iterations import concatenate_log_files, read, plot

if __name__ == "__main__":
    # Path to logfile directory
    path = '/path/to/logfiles'
    # Name of metric to plot. Available: 'f', 'm2', 'rms', 'lambda',
    #  'alpha'
    metric = 'rms'

    # Plotting arguments
    plot_kwargs = {
        # Start of the x-axis. Default is 0.
        'x_start': None,
        # End of the x-axis. Default is len(iterations) - 1.
        'x_end': None,
        # Interval of x-axis ticks. Default is 1.
        'x_interval': None,
        # Start (bottom) of the y-axis. Default is min(values).
        'y_start': None,
        # End of the y-axis. Default is max(values).
        'y_end': None,
        # Interval of y-axis ticks. Default is variance of the values.
        'y_interval': None,
        # Width of the figure in inches. Default is 15.
        'fig_width': None,
        # Height of the figure in inches. Default is 7.5.
        'fig_height': None,
        # Whether or not to add a minor tick mark between each major
        #  tick. Default is False.
        'minor_ticks': False

    }

    logfile = concatenate_log_files(path)
    metrics = read(logfile)
    figure = plot(metric, metrics[metric])
    plotfile = metric + '.png'
    figure.savefig(plotfile)
    print("Complete!")
    print("Concatenated logfile: {}".format(logfile))
    print("Plot: {}".format(plotfile))
