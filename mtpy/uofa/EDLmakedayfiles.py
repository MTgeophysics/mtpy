#!/usr/bin/env python
"""

This is a convenience script for the generation of dayfiles. 
It needs the location of a folder with time series and the sampling period as arguments.

The time series files have to be named in the EDL-ascii output standard, 
which codes stationname and start time of the file in the name. 

The data have to be either single column values or in 2-column form.

        wrapper for the generation of dayfiles for EDL data.

        2 mandatory arguments: 
        - path to files 
        - sampling interval (in seconds)

        3 optional arguments:
        - name of the output directory - cannot start with '-' 
        - stationname - cannot start with '-' 
        - flag '-R (or -r)', if the directory shall be searched for data recursively 

"""


import sys, os
import os.path as op

import mtpy.utils.filehandling as MTfh
reload(MTfh)


def main():

    if len(sys.argv) < 3:
        sys.exit('\nNeed at least 2 arguments: \n\n '
            '<path to files> \n <sampling in seconds> \n\n'
            '[optional: <output dir>] \n [optional: <stationname>]\n'
            '[optional: <recursive flag -R>]\n'
            '(set this option for including all subfolders)\n\n')

    outdir = None
    stationname = None
    recursive = False
    if len(sys.argv) > 3:
        optionals = sys.argv[3:]
        for o in optionals:
            o = o.strip()
            if o[0] == '-':
                if o[1].lower() == 'r':
                    recursive = True
                continue
            elif outdir is None:
                outdir = o
                continue
            elif stationname is None:
                stationname = o 
                continue


    pathname_raw = sys.argv[1]
    pathname = op.abspath(op.realpath(pathname_raw))

    if not op.isdir(pathname):
        sys.exit('Data file(s) path not existing: {0}'.format(pathname))

    try:
        sampling = float(sys.argv[2])
        if sampling <= 0 : raise
    except:
        sys.exit('Second argument must be sampling interval in seconds (int/float)')

    if recursive is True:
        lo_files = []
        for i,j,k in os.walk(pathname):
            lof = [op.abspath(op.join(i,f)) for f in j]
            lo_files.extend(lof)
        pathname = list(set(lo_files))

    MTfh.EDL_make_dayfiles(pathname, sampling, stationname, outdir)



if __name__=='__main__':
    main()
