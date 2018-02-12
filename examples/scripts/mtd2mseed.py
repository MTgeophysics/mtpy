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

Copyright (c) 2014, Lars Krieger/the MTpy team
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""

import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import fnmatch

import mtpy.utils.exceptions as MTex
import mtpy.utils.mseed as MTms
import mtpy.utils.filehandling as MTfh
reload(MTfh)
reload(MTex)
reload(MTms)

def main():

    if len(sys.argv) < 2:
        sys.exit('\n\tNeed at least 1 argument: <path to files> [<output dir>] [<network code>] [<location code>]\n')
        
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
    #we need relative paths here!!!
    indir = pathname_raw   #op.abspath(op.realpath(pathname_raw))

    if not op.isdir(indir):
        raise MTex.MTpyError_inputarguments('Data file(s) path not existing: {0}'.format(indir))

    #define output directory for storing miniSeed files
    #outpath = op.join(os.curdir,'miniSeed')    
    if outdir is not None:
        try:
            outpath = op.abspath(op.join(os.curdir,outdir))
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise
            if not os.access(outpath, os.W_OK):
                raise
        except:
            print 'Cannot generate writable output directory {0} - using generic location "miniSeed" instead'.format(outpath)
            outdir = None
    if outdir is None:
        outpath = op.join(os.curdir,'miniSeed')    
        try: 
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise
            if not os.access(outpath, os.W_OK):
                raise
        except:
            sys.exit('Error ! - Cannot generate writable output directory "miniSeed" - abort...')
    outdir = op.abspath(outpath)
 
    lo_dirs = []
    for i,j,k in os.walk(indir):
        lofolders = [op.join(i,f) for f in j]
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
            outpath = op.abspath(op.join(outdir,i))
            if not op.isdir(outpath):
                os.makedirs(outpath)
            lo_outdirs.append(outpath)
    except:
        raise MTex.MTpyError_inputarguments('ERROR - Cannot set up output directory {0}'.format(outpath))


    for idx_ipath,inpath in enumerate(lo_indirs):
        lo_infiles  = [ i  for i in os.listdir(inpath) if op.isfile(op.abspath(op.join(inpath,i)))]
        lo_outfiles = [ op.abspath(op.join(lo_outdirs[idx_ipath],i)) for i in lo_infiles]
        lo_infiles = [ op.abspath(op.join(inpath,i)) for i in lo_infiles]

       
        for idx_fn, fn in enumerate(lo_infiles):

            if MTfh.validate_ts_file(fn) is False:
                print 'Warning - MT ts data file {0} is not valid (check header)!!!'.format(fn)
                #continue
            print 'reading file {0}'.format(fn)
            outfn = MTms.convertfile_ts2miniseed(fn, lo_outfiles[idx_fn], location=location, network = network)
            print 'wrote file {0}'.format(outfn)

if __name__=='__main__':
    main()
