#!/usr/bin/env python
"""
mtpy/mtpy/uofa/convert_coordinates_in_edis.py

This is a convenience script for converting coordinates in EDI files.

Files are parsed and if a 'lat' or 'lon' is detected, the argument on
the other side of an '=' is converted into decimal degrees. The rest of the file
remains unchanged.

argument:
- directory containing edi files

optional:
- output directory [default: 'decimal_degrees']


"""

import sys
import os
import os.path as op
import fnmatch
import re

import mtpy.utils.format as MTft


def main():

    if len(sys.argv) < 2:
        sys.exit('\nNeed at least 1 arguments:\n '
                 '\n <path to EDI files> \n '
                 '[optional: <output path>]\n')

    edidir = sys.argv[1]
    if not op.isdir(edidir):
        print('Given directory does not exist {0}'.format(edidir))
        sys.exit()

    edilist = []
    try:
        edilist = fnmatch.filter(os.listdir(edidir), '*.[Ee][Dd][Ii]')
        if len(edilist) == 0:
            raise
        edilist = [op.abspath(op.join(edidir, i)) for i in edilist]
    except:
        print('Given directory does not contain edi files: {0}'.format(edidir))

    outputdir = op.join(edidir, 'decimal_degrees')
    if not op.isdir(outputdir):
        os.makedirs(outputdir)

    if len(sys.argv) > 2:
        outputdir = sys.argv[2]
        try:
            if not op.isdir(outputdir):
                os.makedirs(outputdir)
        except:
            print('could not generate output directory - using default')
            outputdir = op.join(edidir, 'decimal_degrees')
            if not op.isdir(outputdir):
                os.makedirs(outputdir)

    path = convert_edis(edilist, outputdir)

    return path


def convert_edis(edilist, output_path):

    for edi in edilist:
        infile = edi
        outfile_raw = os.path.split(edi)[1]
        outfile = op.join(output_path, outfile_raw)
        outstring = ''
        with open(infile, 'r') as F:
            edilines = F.readlines()

        for line in edilines:
            if not ('lat' in line.lower() or 'lon' in line.lower()):
                outstring += line
                continue
            linelist = line.strip().split('=')
            coord = linelist[1]
            dec_coord = str(MTft.assert_decimal_coordinates(coord))

            outstring += '\t{0}={1}\t\n'.format(linelist[0], dec_coord)

        with open(outfile, 'w') as Fout:
            Fout.write(outstring.expandtabs(4))


if __name__ == '__main__':
    main()
