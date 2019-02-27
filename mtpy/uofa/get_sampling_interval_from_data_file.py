#!/usr/bin/env python
"""
    This is a convenience script for finding the sampling interval from a data file.

    It needs the location of a file with time series and optionally the length of the file in seconds; default for the latter is 3600 (hourfile).
"""

import sys
import os.path as op

import mtpy.utils.exceptions as MTex
import mtpy.utils.filehandling as MTfh
# reload(FH)


def main():

    if len(sys.argv) < 2:
        sys.exit('\n usage:\n\t get_sampling_interval_from_data_file.py  <filename>'
                 '[optional: <file length in seconds>]\n\n'
                 'If no second argument is given, a file length of 3600 seconds is assumed.\n\n')

    filename_raw = sys.argv[1]
    filename = op.abspath(op.realpath(filename_raw))

    if not op.isfile(filename):
        raise MTex.MTpyError_inputarguments(
            'File not existing: %s' % (filename))

    try:
        length = float(sys.argv[2])
        if length <= 0:
            raise
    except:
        print('Could not understand second argument - must be a length in seconds (int/float) - set to 3600 ')
        length = 3600

    print(MTfh.get_sampling_interval_fromdatafile(filename, length=length))


if __name__ == '__main__':
    main()
