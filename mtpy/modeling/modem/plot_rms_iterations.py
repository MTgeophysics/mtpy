"""
This module concatenates ModEm log files into a single file and then
plots the RMS across all iterations.

Author: Bren Moushall
"""
import os
import sys
import glob
import re
import argparse

import matplotlib 
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
            metrics[metric].append(float(value))
    return metrics


def plot(metric, values):
    FIGSIZE = 15, 7.5
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_title(metric.upper() + " Across Iterations")
    ax.set_xlabel('Iterations')
    ax.set_xticks(range(1, len(values) - 1, 2))
    ax.set_xticks(range(0, len(values), 1), minor=True)
    ax.set_ylabel(metric.upper())
    interval = np.var(np.asarray(values))
    ax.set_yticks(np.arange(min(values), max(values) + interval, interval))
    ax.set_yticks(np.arange(min(values), max(values) + interval, interval / 2), minor=True)
    ax.set_ylim(min(values), max(values) + interval)
    ax.set_xmargin(0)
    ax.plot(values, color='r', linewidth=2)
    return fig


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """
        Provide a directory containing ModEm logfiles to have them
        joined into a single file. This joint logfile is saved as
        '{directory_name}.log'.

        Logfiles are sorted alphanumerically before they are
        joined, so ensure the files in the directory are named
        correctly to achieve the desired ordering.

        By default, the 'rms' value is plotted, but any value that
        appears in the logfile per iteration as '{metric}={value}'
        can be plotted by provind '-m {metric}'. Plot is saved as
        '{metrc}.png'
        """)
    parser.add_argument('log_directory', help="Path to a directory containing logfiles (can be "
                        "relative)")
    parser.add_argument('--metric', '-m', help="Name of the metric to plot, avaialable are: "
                        "'f', 'm2', 'rms', 'lambda', 'alpha' (default is 'rms')",
                        default='rms')
    args = parser.parse_args()
    logfile = concatenate_log_files(args.log_directory)
    metrics = read(logfile)
    figure = plot(args.metric, metrics[args.metric])
    plotfile = args.metric + '.png'
    figure.savefig(plotfile)
    # TODO: Fix logging
    print("Complete!")
    print("Concatenated logfile: {}".format(logfile))
    print("Plot: {}".format(plotfile))
