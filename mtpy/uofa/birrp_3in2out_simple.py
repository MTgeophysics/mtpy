#!/usr/bin/env python
"""
This is a convenience script for running BIRRP.

arguments:
birrp executable, stationname (uppercase), directory containing time series files, coherence threshold

A subfolder 'birrp_processed' for the output is generated within the time series directory

"""


import numpy as np
import re
import sys
import os
import glob
import os.path as op
import glob
import calendar
import time


import mtpy.utils.exceptions as MTex

import mtpy.processing.birrp as MTbp
#reload(MTbp)


def main():

    if len(sys.argv) < 4:
        print('\nNeed at least 3 arguments: <path to BIRRP executable> '\
            '<station name> <directory for time series>\n\n'\
            'Optional arguments: \n [coherence threshold]\n'\
            ' [start time] \n [end time]\n\n')
        return

    try:
        coherence_th = float(sys.argv[4])
        if not 0 <= coherence_th <= 1:
            raise
    except:
        print('coherence value invalid (float from interval [0,1]) - set to 0 instead')
        coherence_th = 0

    try:
        starttime = float(sys.argv[5])
    except:
        starttime = None

    try:
        endtime = float(sys.argv[6])
    except:
        endtime = None

    birrp_exe_raw = sys.argv[1]
    birrp_exe = op.abspath(op.realpath(birrp_exe_raw))

    if not op.isfile(birrp_exe):
        raise MTex.MTpyError_inputarguments(
            'Birrp executable not existing: %s' %
            (birrp_exe))

    stationname = sys.argv[2].upper()

    ts_dir_raw = sys.argv[3]
    ts_dir = op.abspath(op.realpath(ts_dir_raw))

    if not op.isdir(ts_dir):
        raise MTex.MTpyError_inputarguments(
            'Time series directory not existing: %s' %
            (ts_dir))

    if 1:
        MTbp.runbirrp_Nin2out_simple(birrp_exe, stationname, ts_dir, coherence_th,
                                     None, 3, None, starttime, endtime)
    # except:
    #     print 'ERROR - Could not process input data using BIRRP'


if __name__ == '__main__':
    main()
