#!/usr/bin/env python
"""

This is a convenience script for the generation of dayfiles.
It needs the location of a folder with time series and the sampling period as arguments.

The time series files have to be named in the EDL-ascii output standard,
which codes stationname and start time  of the file in the name.

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


import sys
import os
import os.path as op
import mtpy.utils.exceptions as MTex

import mtpy.utils.filehandling as MTfh
#reload(MTfh)


def main():

    if len(sys.argv) < 3:
        sys.exit('\nNeed at least 4 arguments: \n\n '
                 '<path to files> \n <sampling in seconds> \n'
                 '<output dir> \n <stationname>\n'
                 '[optional: <recursive flag -R>]\n'
                 '(set this option for including all subfolders)\n\n')

    print()

    outdir = None
    stationname = None
    recursive = False

    multiple_stations = False

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

    if stationname is not None:
        # check, if it's actually a comma-separated list:
        try:
            stationlist = stationname.split(',')
            if len(stationlist) > 1:
                multiple_stations = True
                stationlist = [i.upper() for i in stationlist]
        except:
            stationlist = [stationname]
    else:
        stationlist = [None]

    print(stationlist)

    pathname_raw = sys.argv[1]
    pathname = op.abspath(op.realpath(pathname_raw))

    if not op.isdir(pathname):
        sys.exit('Data file(s) path not existing: {0}\n'.format(pathname))

    try:
        sampling = float(sys.argv[2])
        if sampling <= 0:
            raise
    except:
        sys.exit('Second argument must be sampling interval in seconds (int/float)')

    if recursive is True:
        lo_folders = []
        for i, j, k in os.walk(pathname):
            lof = [op.abspath(op.join(i, f)) for f in j]
            if stationname is not None:
                for stationname in stationlist:
                    for curr_folder in lof:
                        content_of_folder = os.listdir(curr_folder)
                        # print curr_folder
                        lof_station = [
                            i for i in content_of_folder if stationname.lower() in i.lower()]
                        if len(lof_station) > 0:
                            lo_folders.append(curr_folder)
        pathname = list(set(lo_folders))

    if len(pathname) == 0:
        sys.exit(
            '\n\tERROR - No (sub-) folders for stations {0} found\n'.format(stationlist))

    for stationname in stationlist:
        print('....\n')
        print('processing station ', stationname.upper())
        # if pathname[0] is not None:
        #     station_pathname = [i for i in pathname if stationname.lower() in i.lower()]
        #     if len(station_pathname) == 0:
        #         station_pathname = None
        # else:
        station_pathname = pathname

        try:
            MTfh.EDL_make_dayfiles(
                station_pathname,
                sampling,
                stationname.upper(),
                outdir)
        except MTex.MTpyError_inputarguments:
            if stationname is None:
                sys.exit('\n\tERROR - No data found in (sub-)folders\n')
            else:
                sys.exit(
                    '\n\tERROR - No data found in (sub-)folders for station {0}\n'.format(stationname.upper()))
        except MemoryError:
            sys.exit('\n\tERROR - Not enough memory to store temporary arrays!\n')
        except:
            sys.exit('\n\tERROR - could not process (sub-)folders')

    print('\n')


if __name__ == '__main__':
    main()
