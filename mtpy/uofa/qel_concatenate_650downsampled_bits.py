#!/usr/bin/env python


import os
import sys
import os.path as op

import numpy as np


def main():

    if len(sys.argv) < 3:
        sys.exit('\n\tERROR - need 2 arguments as input: '
                 '<input files location> <output file name> \n')

    indir = sys.argv[1]
    outfilename = sys.argv[2]

    indir = op.join(op.abspath(os.curdir), indir)
    print()
    if not op.isdir(indir):
        print('WARNING - no such directory: %s\n' % (indir))
        sys.exit()

    outfn = op.join(op.abspath(os.curdir), outfilename)

    concatenate1station650(indir, outfn)

    print('Downsampled data written to file %s ' % (outfn))

    print()


def concatenate1station650(indir, outfile):
    """
    Concatenate ascii file entries into a single 5-column first differences file.

    The header line contains starting values, all later lines show differences
    to the relative preceding line.

    """

    lo_files = sorted(os.listdir(indir))
    lo_files = [op.join(indir, i) for i in lo_files]

    if len(lo_files) == 0:
        print('ERROR - no file in directory %s\n' % (indir))
        sys.exit()

    # open output file
    Fout = open(outfile, 'w')

    # initialise carrying varibles
    header = False
    N = None
    E = None
    S = None
    W = None
    t0 = None

    print('Browsing/reading files in %s' % (indir))

    # go through files in directory
    for infile in lo_files:
        try:
            Fin = open(infile)
            # scan all lines
            for line in Fin:

                indata = line.strip().split()
                # 4 channel values are integers
                indata[:-1] = [int(float(i)) for i in indata[:-1]]
                # 1 time stamp is high precision float64
                indata[-1] = np.round(np.float64(indata[-1]), 2)

                if header is False:
                    # set first values of first file as initial points
                    N = int(indata[1])
                    E = int(indata[2])
                    S = int(indata[3])
                    W = int(indata[4])
                    t0 = int(indata[5])

                    headerline = '# %d\t%d\t%d\t%d\t%d\n' % (N, E, S, W, t0)

                    # write to header line
                    Fout.write(headerline)

                    # one header per file:
                    header = True

                # define ongoing line entries as differences to values from
                # line before
                currentline = '%d\t%d\t%d\t%d\t\n' % (
                    indata[1] - N, indata[2] - E, indata[3] - S, indata[4] - W)
                Fout.write(currentline)

                # update carrying variables
                N = int(indata[1])
                E = int(indata[2])
                S = int(indata[3])
                W = int(indata[4])

            Fin.close()

        except:
            print('\tWARNING - could not read data from file: %s' % (infile))
            continue

    Fout.close()


if __name__ == '__main__':
    main()
