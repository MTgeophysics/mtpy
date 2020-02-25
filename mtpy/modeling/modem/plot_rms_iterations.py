"""
This module concatenates ModEm log files into a single file and then
plots the RMS across all iterations.

Author: Bren Moushall
"""
import os
import sys
import glob
import re

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
            all_lines.extend(f.readlines())

    new_logfile = os.path.basename(directory) + '.log'
    new_logfile = os.path.join(directory, new_logfile)
    with open(new_logfile, 'w') as f:
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
    value_lines = [l.strip().replace('   ', '').split() for l in value_lines]
    print(len(value_lines))
    metrics = {'f': [], 'm2': [], 'rms': [], 'lambda': [], 'alpha': []}
    for line in value_lines:
        del line[0]  # Get rid of 'with:'
        for word in line:
            metric, value = word.split('=')
            metrics[metric].append(value)
    return metrics


def plot(metric, values):
    FIGSIZE = 25, 12.5
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlabel('Iterations')
    ax.set_ylabel(metric)
    values = list(reversed(values))
    ax.plot(np.arange(1, len(values) + 1), values)
    return fig


if __name__ == '__main__':
    logfile = concatenate_log_files(sys.argv[1])
    metrics = read(logfile)
    figure = plot(sys.argv[2], metrics[sys.argv[2]])
    figure.savefig(sys.argv[2] + '.png')

