#! /usr/bin/env python
"""
Description:
    Utility functions for unit-testing suite


CreationDate:   2/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate: 07/11/2017  YG
                02/11/2017  FZ
"""

# package tests.beta scope global params
import os
import shutil
import sys
from difflib import unified_diff

import matplotlib

from tests import TEST_MTPY_ROOT

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt

    plt.ion()

EDI_DATA_DIR = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files'))
EDI_DATA_DIR2 = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files_2'))

AUS_TOPO_FILE = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/AussieContinent_etopo1.asc'))

# path to directory containing model input files - samples reference for compare
SAMPLE_DIR = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/model_files'))  # r'E:\Githubz\mtpy\examples\model_files'


def _diff_files(after, before):
    """
    compare two files using diff
    :param before:
    :param after:
    :return: the number count of different lines
    """

    with open(before) as f2p:
        before_lines = f2p.readlines()
    with open(after) as f1p:
        after_lines = f1p.readlines()

    msg = "Comparing {} and {}:\n".format(before, after)
    is_identical = False

    lines = [line for line in unified_diff(
        before_lines,
        after_lines,
        fromfile="baseline ({})".format(before),
        tofile="test ({})".format(after),
        n=0)]
    if lines:
        msg += "  Found differences:\n    " + "    ".join(lines)
        is_identical = False
    else:
        msg += " NO differences found."
        is_identical = True

    return is_identical, msg



def _clean_recreate(adir):
    if os.path.exists(adir):
        # clear dir if it already exist
        shutil.rmtree(adir)

    os.mkdir(adir)
