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
import os

from mtpy.modeling.modem.plot_rms_iterations import concatenate_log_files, read, plot

if __name__ == "__main__":
    # Path to logfile directory
    path = '/home/bren/data_mtpy/rms_logs'
    # Name of metric to plot. Available: 'f', 'm2', 'rms', 'lambda',
    #  'alpha'
    metric = 'rms'
    # Plotting arguments
    plot_kwargs = {
        # Set as None to use default values.
        # Interval of x-axis ticks - default is 1
        'x_interval': 5.,
        # Interval of y-axis ticks - default is variance of the values
        'y_interval': 0.5,
        # Width of the figure in pixels - default is 800
        'fig_width': 800,
        # Height of the figure in inches - default is 800
        'fig_height': 800,
        # DPI of the figure - default is 100
        'dpi': None
    }

    logfile = concatenate_log_files(path)
    metrics = read(logfile)
    figure = plot(metric, metrics[metric], **plot_kwargs)
    plotfile = os.path.join(path, metric + '.png')
    figure.savefig(plotfile)
    print("Complete!")
    print("Concatenated logfile: {}".format(logfile))
    print("Plot: {}".format(plotfile))
