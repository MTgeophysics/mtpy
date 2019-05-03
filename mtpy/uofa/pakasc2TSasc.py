#!/usr/bin/env python


import sys
import os
import os.path as op
import numpy as np
import glob


def pakascii2TSascii(fn):
    """
    ...

    """
    try:
        F = open(fn, 'r')
        raw_data = F.readlines()
        F.close()
    except:
        sys.exit('ERROR - input file could not be read:\n%s' % (fn))

    lo_datalines = []
    n_samples = len(raw_data) - 4

    if n_samples < 1:
        sys.exit('ERROR - no data in PakAscii file:\n%s' % (fn))

    if 1:
        for step in range(n_samples):
            current_line_idx = 4 + step
            current_line = raw_data[current_line_idx].strip().split()
            if raw_data[current_line_idx].strip()[0] == '#':
                continue
            ch0 = float(current_line[1])
            ch1 = float(current_line[2])
            ch2 = float(current_line[3])
            ch3 = float(current_line[4])
            t = float(current_line[12])
            lo_datalines.append([t, ch0, ch1, ch2, ch3])
    else:
        sys.exit('ERROR -  PakAsciifile contains errorneous line')

    data_array = np.array(lo_datalines)

    outfn = op.splitext(op.basename(fn))[0] + '.dat'

    np.savetxt(outfn, data_array)

    return outfn


def main():

    arglist = sys.argv[1:]

    filelist_raw = []

    for arg in arglist:
        globlist = glob.glob(arg)
        if len(globlist) == 0:
            print('Warning -- cannot read given file(s):\n  %s' % (arg))
        filelist_raw.extend(globlist)

    filelist = filelist_raw

    # for fr in filelist_raw:
    #     try:
    #         if op.isfile(op.abspath(fr)):
    #             filelist.append(fr)
    #         else:
    #             raise
    #     except:
    #         print 'Warning -- cannot read given file:\n  %s'%(fr)
    #         continue

    if len(filelist) == 0:
        sys.exit('ERROR - no files given for conversion')

    for f in filelist:
        try:
            print('converting file %s ...' % (f))
            pakascii2TSascii(f)

        except:
            print('Warning -- cannot convert given file:\n  %s' % (f))
            continue

    return 0


#===================================================

if __name__ == "__main__":
    main()
