#!/usr/bin/env python
"""
    This is a convenience script for finding the sampling interval from a data file. 

    It needs the location of a file with time series and optionally the length of the file in seconds; default for the latter is 3600 (hourfile).
"""

from mtpy.utils.exceptions import *
import sys
import os.path as op

import mtpy.utils.filehandling as FH
reload(FH)



def main():

    if len(sys.argv) < 2:
        raise MTpyError_inputarguments('Need at least 1 argument: <filename>')


    filename_raw = sys.argv[1] 
    filename = op.abspath(op.realpath(filename_raw))


    if not op.isfile(filename):
        raise MTpyError_inputarguments('File not existing: %s' % (filename))

    try:
        length = float(sys.argv[1])
        if length <= 0 : raise
    except:
        print 'Could not understand second argument - must be a length in seconds (int/float) - set to 3600 '
        length = 3600

    print FH.get_sampling_interval_fromdatafile(filename, length = length)


if __name__=='__main__':
    main()
