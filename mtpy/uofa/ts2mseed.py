#!/usr/bin/env python
"""

    This is a convenience script for the conversion of MTpy-style time series data files to miniSeed.

    It needs the location of a folder with TS data. The folder structure is mimicked and the TS converted into miniSeed. The folder structure is then put into the destination directory.
    The latter can be given as an argument, if not, a local directory 'miniSeed' will be generated in the current directory


        1 mandatory argument:
        - path to TS data files


        3 optional arguments (in this order):
        - name of the output directory
        - 'location' code (2 characters)
        - 'network' code (2 characters)

"""

import numpy as np
import re
import sys
import os
import glob
import os.path as op
import glob
import fnmatch

import mtpy.utils.exceptions as MTex
import mtpy.utils.mseed as MTms
import mtpy.utils.filehandling as MTfh
#reload(MTfh)
#reload(MTex)
#reload(MTms)


def main():

    if len(sys.argv) < 2:
        sys.exit(
            '\n\tNeed at least 1 argument: <path to files> [<output dir>] [<network code>] [<location code>]\n')

    outdir = None
    location = ''
    network = ''

    if len(sys.argv) > 2:
        outdir = sys.argv[2]
        if len(sys.argv) > 3:
            network = sys.argv[3]
            if len(sys.argv) > 4:
                location = sys.argv[4]

    pathname_raw = sys.argv[1]
    # we need relative paths here!!!
    indir = pathname_raw  # op.abspath(op.realpath(pathname_raw))

    if not op.isdir(indir):
        raise MTex.MTpyError_inputarguments(
            'Data file(s) path not existing: {0}'.format(indir))

    # define output directory for storing miniSeed files
    #outpath = op.join(os.curdir,'miniSeed')
    if outdir is not None:
        try:
            outpath = op.abspath(op.join(os.curdir, outdir))
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise
            if not os.access(outpath, os.W_OK):
                raise
        except:
            print('Cannot generate writable output directory {0} - using generic location "miniSeed" instead'.format(outpath))
            outdir = None
    if outdir is None:
        outpath = op.join(os.curdir, 'miniSeed')
        try:
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise
            if not os.access(outpath, os.W_OK):
                raise
        except:
            sys.exit(
                'Error ! - Cannot generate writable output directory "miniSeed" - abort...')
    outdir = op.abspath(outpath)

    lo_dirs = []
    for i, j, k in os.walk(indir):
        lofolders = [op.join(i, f) for f in j]
        lo_dirs.extend(lofolders)
    lo_dirs.append(indir)
    pathname = list(set(lo_dirs))
    if len(pathname) == 0:
        pathname = [indir]

    lo_indirs = pathname
    lo_outdirs = []
    #'pathname' is a list of relative pathnames. to be reconstructed under the given 'outdir'
    try:
        for i in lo_indirs:
            outpath = op.abspath(op.join(outdir, i))
            if not op.isdir(outpath):
                os.makedirs(outpath)
            lo_outdirs.append(outpath)
    except:
        raise MTex.MTpyError_inputarguments(
            'ERROR - Cannot set up output directory {0}'.format(outpath))

    for idx_ipath, inpath in enumerate(lo_indirs):
        lo_infiles = [
            i for i in os.listdir(inpath) if op.isfile(
                op.abspath(
                    op.join(
                        inpath,
                        i)))]
        lo_outfiles = [
            op.abspath(
                op.join(
                    lo_outdirs[idx_ipath],
                    i)) for i in lo_infiles]
        lo_infiles = [op.abspath(op.join(inpath, i)) for i in lo_infiles]

        for idx_fn, fn in enumerate(lo_infiles):

            if MTfh.validate_ts_file(fn) is False:
                print('Warning - MT ts data file {0} is not valid (check header)!!!'.format(fn))
                # continue
            print('reading file {0}'.format(fn))
            outfn = MTms.convertfile_ts2miniseed(
                fn, lo_outfiles[idx_fn], location=location, network=network)
            print('wrote file {0}'.format(outfn))

if __name__ == '__main__':
    main()
