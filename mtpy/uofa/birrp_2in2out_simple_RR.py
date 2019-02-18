#!/usr/bin/env python
"""
Convenience script for running BIRRP with one remote reference.

arguments:
birrp executable, stationname, remote refernce station name, directory containing time series files

optional:
coherence threshold, start time of processing, end time of processing

A subfolder 'birrp_processed_rr' for the output is generated within the time series directory

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

    if len(sys.argv) < 5:
        print('\nNeed at least 4 arguments: <path to BIRRP executable> '\
            '<station name> <remote station name> <directory for'\
            ' time series>\n\n'\
            'Optional arguments: \n [coherence threshold]\n'\
            ' [start time] \n [end time]\n\n')
        return

    try:
        coherence_th = float(sys.argv[5])
        if not 0 < coherence_th <= 1:
            raise
    except:
        print(' Warning - Coherence value invalid (float from interval ]0,1]) '\
            '- set to 0.5 instead')
        coherence_th = 0.5

    try:
        starttime = float(sys.argv[6])
    except:
        starttime = None

    try:
        endtime = float(sys.argv[7])
    except:
        endtime = None

    birrp_exe_raw = sys.argv[1]
    birrp_exe = op.abspath(op.realpath(birrp_exe_raw))

    if not op.isfile(birrp_exe):
        print('\nError - Birrp executable not existing: {0}\n'.format(birrp_exe))

    stationname = sys.argv[2].upper()
    rr_stationname = sys.argv[3].upper()

    ts_dir_raw = sys.argv[4]
    ts_dir = op.abspath(op.realpath(ts_dir_raw))

    if not op.isdir(ts_dir):
        print('\nError - Time series directory not existing: {0}\n'.format(ts_dir))
        return

    if 1:
        MTbp.runbirrp2in2out_simple(birrp_exe, stationname, ts_dir, coherence_th,
                                    rr_stationname, None, starttime, endtime)
    # except:
    #     print '\n\tERROR - Could not process input data using BIRRP\n'
    #     return


if __name__ == '__main__':
    main()
