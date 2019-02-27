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
        sys.exit('\n\tNeed at least 1 argument:\n\n <path to files>\n[optional:'
                 '<output dir>] \n[optional: <combine channels option -c>]\n'
                 '[optional: <invert West and South option -i>]\n\n')

    outdir = None
    combine_flag = False
    inversion_flag = False

    if len(sys.argv) > 2:
        optionals = sys.argv[2:]

        for idx_o, o in enumerate(optionals):
            if o[0] != '-':
                outdir = o
                continue
            option = o[1].lower()
            if option not in ['c', 'i']:
                print('unknown option: {0}'.format(option))
                continue
            if option == 'c':
                combine_flag = True
            if option == 'i':
                inversion_flag = True

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
            print('Cannot generate writable output directory {0} - using generic'\
                ' location "ascii" instead'.format(outpath))
            outdir = None
    if outdir is None:
        outpath = op.join(os.curdir, 'ascii')
        try:
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise
            if not os.access(outpath, os.W_OK):
                raise
        except:
            sys.exit('Error ! - Cannot generate writable output directory '
                     '"ascii" - abort...')
    outdir = op.abspath(outpath)

    convert_ms2ts(indir, outdir, combine_flag, inversion_flag)


def convert_ms2ts(indir, outdir, combine=True, invert=False):

    lo_dirs = []
    for i, j, k in os.walk(indir):
        lofolders = [op.join(i, f) for f in j]
        lo_dirs.extend(lofolders)
    lo_dirs.append(indir)
    pathname = list(set(lo_dirs))
    if len(pathname) == 0:
        pathname = [indir]

    lo_indirs = [op.relpath(i, start=indir) for i in pathname]
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
    lo_indirs = [op.join(indir, i) for i in lo_indirs]

    for idx_ipath, inpath in enumerate(lo_indirs):
        lo_infiles = [i for i in os.listdir(inpath) if
                      op.isfile(op.abspath(op.join(inpath, i)))]

        lo_outfiles = [op.abspath(op.join(lo_outdirs[idx_ipath], i)) for
                       i in lo_infiles]

        lo_infiles = [op.abspath(op.join(indir, inpath, i))
                      for i in lo_infiles]

        for idx_fn, fn in enumerate(lo_infiles):

            print('reading file {0}'.format(fn))
            try:
                outfn = MTms.quadrupol_convertfile_miniseed2ts(fn,
                                                               lo_outfiles[idx_fn], combine, invert)

                print('wrote file(s) {0}'.format(outfn))
            except:
                print('Warning - file {0} is not in valid miniseed  format!!!'.format(fn))
                continue


if __name__ == '__main__':
    main()
