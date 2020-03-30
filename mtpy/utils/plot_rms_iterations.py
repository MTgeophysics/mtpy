from __future__ import division
import os
import glob
import re

import matplotlib.pyplot as plt
import numpy as np


def _an_sort(collection):
    """Alphanumeric sort.
    """
    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(collection, key=alphanum_key)


def concatenate_log_files(directory):
    """
    Any file of the pattern '*.log' will be included. The files are
    sorted alphanumerically and this is the order they will be
    concatenated in. It is up to the user to ensure that the files are
    named correctly to achieve the desired order.

    Args:
        directory (str):
    """
    directory = os.path.abspath(directory)
    files = _an_sort(glob.glob(os.path.join(directory, '*.log')))
    if not files:
        raise ValueError("No log files found in '{}'.".format(directory))
    all_lines = []
    for fn in files:
        with open(fn, 'r') as f:
            lines = f.readlines()
            # Skip already joined logfiles
            if lines[0] == 'Concatenated log files\n':
                continue
            else:
                all_lines.extend(lines)

    new_logfile = os.path.basename(directory) + '.log'
    new_logfile = os.path.join(directory, new_logfile)
    all_lines.insert(0, 'Concatenated log files\n')
    with open(new_logfile, 'w+') as f:
        f.writelines(all_lines)
    return new_logfile


def read(logfile):
    """
    Get a sequence of values from a ModEM logfile. Each type of value
    present in the logfile is collected and ordered by iteration.


    Args:
        path (str): Path to the logfile to be read.

    Returns
        dict of str, float: A dictionary containing lists of metric
            values.
    """
    lf = open(logfile, 'r')
    lines = lf.readlines()
    value_lines = [l for l in lines if l.strip().startswith('with:')]
    if not value_lines:
        raise ValueError("Concatenated log file did not contain any readable values")
    value_lines = [l.strip().replace('   ', '').split() for l in value_lines]
    metrics = {'f': [], 'm2': [], 'rms': [], 'lambda': [], 'alpha': []}
    for line in value_lines:
        del line[0]  # Get rid of 'with:'
        for word in line:
            metric, value = word.split('=')
            try:
                value = float(value)
            except ValueError:
                value = None
            metrics[metric].append(value)
    return metrics


def plot(metric, values, x_start=0, x_end=None, x_interval=1, y_start=None, y_end=None,
         y_interval=None, fig_width=1900, fig_height=1200, dpi=100, minor_ticks=True):
    fig_width = 800 if fig_width is None else fig_width
    fig_height = 800 if fig_height is None else fig_height
    dpi = 100 if dpi is None else dpi
    # Convert pixels to inches
    figsize = 0.0104166667 * fig_width, 0.0104166667 * fig_height
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_title(metric.upper() + "/Iteration")

    ax.set_xlabel('Iterations')
    x_start = 0 if x_start is None else x_start
    x_end = len(values) if x_end is None else x_end
    x_interval = 1 if x_interval is None else x_interval
    ax.set_xticks(np.arange(x_start, x_end, x_interval))
    if minor_ticks:
        ax.set_xticks(np.arange(x_start, x_end, x_interval / 2), minor=True)

    ax.set_ylabel(metric.upper())
    y_start = min(values) if y_start is None else y_start
    y_end = max(values) if y_end is None else y_end
    y_interval = np.var(np.asarray(values)) if y_interval is None else y_interval
    ax.set_yticks(np.arange(1, y_end, y_interval))
    if minor_ticks:
        ax.set_yticks(np.arange(1, y_end, y_interval / 2), minor=True)

    ax.plot(values, color='r', linewidth=2)

    # Set y-lim based on user limit and default padding
    y_lim_bottom = 1 - abs(min(values) - ax.get_ylim()[0])
    y_lim_top = y_end + abs(max(values) - ax.get_ylim()[1])
    ax.set_ylim(y_lim_bottom, y_lim_top)
    return fig



