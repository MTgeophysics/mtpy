#! /usr/bin/env python
"""
Description:
    Utility functions for unit-testing suite


CreationDate:   2/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     2/11/2017   FZ
"""

# package tests.beta scope global params
import os
import shutil

import matplotlib.pyplot as plt
plt.ion()  # all beta/test_ plot disappear automatically

# MTPY_ROOT='/Softlab/Githubz/mtpy'    # source code root dir
MTPY_ROOT='E:/Githubz/mtpy'    # source code root dir

EDI_DATA_DIR = os.path.join(MTPY_ROOT,'examples/data/edi_files')
EDI_DATA_DIR2 = os.path.join(MTPY_ROOT,'examples/data/edi_files_2')

AUS_TOPO_FILE = os.path.join(MTPY_ROOT,'examples/data/AussieContinent_etopo1.asc')

# path to directory containing model input files - samples reference for compare
SAMPLE_DIR = os.path.join(MTPY_ROOT,'examples/model_files') # r'E:\Githubz\mtpy\examples\model_files'

# test runs output directory where to save plots/files to
TEMP_OUT_DIR = os.path.join(MTPY_ROOT,'temp/beta_out_dir')  # r'E:\Githubz\mtpy\temp'


def _diffiles(f1, f2):
    """
    compare two files line-by-line
    :param f1:
    :param f2:
    :return: the number count of different lines
    """

    count = 0

    with open(f1) as f1p:
        test_lines = f1p.readlines()
    with open(f2) as f2p:
        correct_lines = f2p.readlines()

    for test, correct in zip(test_lines, correct_lines):
        if test != correct:
            # print ("Diffiles() Failure@: Expected %r; BUT Got %r." % (correct, test))
            count = count + 1
        else:
            pass

    return count


def _clean_recreate(adir):
    if os.path.exists(adir):
        # clear dir if it already exist
        shutil.rmtree(adir)

    os.mkdir(adir)
