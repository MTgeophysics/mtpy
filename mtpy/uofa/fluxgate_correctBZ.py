#!/usr/bin/env python

"""Correcting the instrumentally hardcoded damping factor 2.2 from the data.

Applies only to the Z component of long period Fluxgate magnetic sensors!!!

Works in-place - so make sure, it's only run once!!!!

"""


import os
import sys
import os.path as op
import numpy as np

import ipdb

if len(sys.argv) < 2:
    sys.exit('\nNeed at least 1 argument: \n\n '
             '<path to files> \n \n'
             '[optional: <output dir>]')


outdir = None
filepath = sys.argv[1]

if len(sys.argv) > 2:
    outpath = sys.argv[2]

pathname = op.abspath(op.realpath(filepath))

if not op.isdir(pathname):
    sys.exit('Data file(s) path not existing: {0}'.format(pathname))

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
    print('Cannot generate writable output directory {0} - using'\
        ' generic location "dayfiles" instead'.format(outdir))
    outpath = pathname
    pass


lo_files = os.listdir(pathname)

lo_files = [i for i in lo_files if i.lower().endswith('.bz')]

if len(lo_files) == 0:
    sys.exit('ERROR - no BZ data in directory {0} \n'.format(pathname))


print()


for bz_file in lo_files:

    try:
        infile = op.join(pathname, bz_file)

        with open(infile) as F:
            header = F.readline()

        data = np.loadtxt(infile)

        outfile = os.path.join(outpath, bz_file)

        Fout = open(outfile, 'w')
        Fout.write(header)
        np.savetxt(Fout, 2.2 * data)
        Fout.close()

        print('\tcorrecting file {0} ....output: {1}'.format(infile, outfile))

    except:
        print('\n\t\tERROR - could not correct file {0}\n'.format(bz_file))

print('\n...Done\n')
