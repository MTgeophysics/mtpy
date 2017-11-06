#! /usr/bin/env python
"""
Description:
    Utility functions for unit-testing suite


CreationDate:   2/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     2/11/2017   FZ
"""

# import section
import os
import shutil

def diffiles(f1, f2):
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

def clean_recreate(adir):
    if os.path.exists(adir):
        # clear dir if it already exist
        shutil.rmtree(adir)

    os.mkdir(adir)
