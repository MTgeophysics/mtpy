#!/usr/bin/env python
"""This is a convenience script for the generation of dayfiles. 
It needs the location of a folder with time series and the sampling period as arguments.

The time series files have to be named in the EDL-ascii output standard, 
which codes stationname and start time of the file in the name. 

The data have to be either single column values or in 2-column form.

"""






import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time


from mtpy.utils.exceptions import *


import mtpy.processing.filehandling as FH
reload(FH)







def main():

    if len(sys.argv) < 3:
        raise MTpyError_inputarguments('Need 2 arguments: <path to files> <sampling in seconds>')


    pathname_raw = sys.argv[1] 
    pathname = op.abspath(op.realpath(pathname_raw))


    if not op.isdir(pathname):
        raise MTpyError_inputarguments('Path not existing: %s' % (pathname))

    try:
        sampling = float(sys.argv[2])
        if sampling <= 0 : raise
    except:
        raise MTpyError_float('Second argument must be sampling interval in seconds (int/float)')


    FH.EDL_make_dayfiles(pathname, sampling)



if __name__=='__main__':
    main()
