#!/usr/bin/env python

"""
Convert an EDI file to multi-column ascii file

first column is frequency, the rest goes:

XXreal, XXimag XXsigma, XYreal, XYimag, XYsigma, ...


"""


import os
import sys
import os.path as op
import numpy as np

import mtpy.core.edi as MTedi


def main():
    print()
    if len(sys.argv) < 2:
        sys.exit('\nNeed at least 1 argument:\n '
                 '<EDI file> \n\n')

    edifn = sys.argv[1]

    try:
        edifn = op.realpath(op.abspath(op.join(os.curdir, edifn)))
        if not op.isfile(edifn):
            raise
    except:
        print('\tERROR - EDI file does not exist: {0}\n'.format(edifn))
        sys.exit()

    try:
        convert2columns(edifn)
    except:
        print('\tERROR - could not convert EDI file\n')


def convert2columns(fn):

    try:
        e_object = MTedi.Edi(filename=fn)
    except:
        print('invalid EDI file')
        raise

    filebase = op.splitext(fn)[0]
    outfn = filebase + '.columns'
    outstring = ''
    for i, f in enumerate(e_object.freq):
        outstring += '{0:> 10.5f}  '.format(f)
        for j in range(2):
            for k in range(2):
                outstring += '\t{0:.4f}  {1:.4f}  {2:.4f}  '.format(float(np.real(e_object.Z.z[i, j, k])),
                                                                    float(np.imag(e_object.Z.z[i, j, k])), e_object.Z.z_err[i, j, k])
        outstring += '\n'

    with open(outfn, 'w') as F:
        F.write(outstring)

    print('\tWritten data to file: {0}\n'.format(outfn))

    return outfn


if __name__ == '__main__':
    main()
