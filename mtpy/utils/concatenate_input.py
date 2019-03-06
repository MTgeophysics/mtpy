#!/usr/bin/env python
"""
Description:
    This script collates data from raw data files in a folder, within a time-range
    provided by the user and outputs corresponding .EX, .EY, .EZ, .BX, .BY and .BZ
    files in an output folder.

References:

CreationDate:   2017/10/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     2017/10/23   RH
"""

import os
import re
import glob
from collections import defaultdict
import numpy as np
import click
from datetime import datetime

class Data():
    def __init__(self, dataPath, startDateTime='', endDateTime=''):
        """

        :param dataPath: path to input files
        :param startDateTime: start date and time
        :param endDateTime: end data and time
        """

        def createPythonDateTime(strDateTime):
            """
            Converts a date-time in string format to a python datetime object
            :param strDateTime: date-time in string format
            :return: datetime object
            """

            if(strDateTime==''): return None

            year = 0
            month = 1
            day = 1
            hour = 0
            minute = 0
            second = 0

            if(len(strDateTime) >= 4): year = int(strDateTime[0:4])
            if(len(strDateTime) >= 6): month = int(strDateTime[4:6])
            if(len(strDateTime) >= 8): day = int(strDateTime[6:8])
            if(len(strDateTime) >= 10): hour = int(strDateTime[8:10])
            if(len(strDateTime) >= 12): minute = int(strDateTime[10:12])
            if(len(strDateTime) >= 14): second = int(strDateTime[12:14])

            dt = datetime(year, month, day, hour, minute, second)
            return dt
        #end func

        # Collect list of files
        self._files = glob.glob(os.path.join(dataPath, '*.*'))

        assert len(self._files), 'Error: No files found. Aborting..'

        # Create date-time objects for range provided
        self._startDateTime = createPythonDateTime(startDateTime)
        self._endDateTime = createPythonDateTime(endDateTime)

        self._filesDict = defaultdict(str)
        for f in self._files:
            dateString = re.findall(r'\d+', os.path.basename(f))

            assert len(dateString) == 1, 'Error processing file name %s. Aborting..'%(f)
            dt = createPythonDateTime(dateString[0])

            #Filter files found based on user-defined time-range
            if(self._startDateTime is not None and dt < self._startDateTime):
                continue
            # end if
            if (self._endDateTime is not None and dt > self._endDateTime):
                continue
            # end if

            self._filesDict[dt] = f
            print(('Selecting %s ..'%os.path.basename(f)))
        #end for
        print(('\nAdded %d files'%(len(list(self._filesDict.keys())))))
    #end func

    def ouput(self, prefix, outputPath):
        """
        :param prefix: output file prefix
        :param outputPath: output folder
        """

        def readData(fileName):
            """
            This function needs to be capable of handling ASCII data files in various formats.
            So far, it can deal with the following variants:
                1. File with no headers and 24 columns
                2: File with 13 lines of header and 7 columns

            :param fileName: Name of the input file
            :return: 3 column b and e fields
            """

            def isHeader(line):
                ssl = re.findall(r'[a-zA-Z]+', line)
                for ss in ssl:
                    # Data lines may have single characters e.g. E, S, etc. for orientation.
                    # Consider a line to be header if there are character strings of length 2 and above
                    if(len(ss)>2): return 1
                # end for

                return 0
            # end func

            def numHeaderLines(fileName):
                f = open(fileName)
                result = 0
                for line in f:
                    h = isHeader(line)
                    if(h == 0): break
                    result += h
                # end for
                f.close()

                return result
            #end func

            # Check for headers
            nhl = numHeaderLines(fileName)

            b = None
            e = None
            if(nhl==0):
                # File variant 1

                d = np.genfromtxt(fileName, invalid_raise=False)
                b = d[:, 6:9]
                e = d[:, 11:14]
            elif(nhl==13):
                # File variant 2

                d = np.genfromtxt(fileName, invalid_raise=False, skip_header=nhl)
                b = d[:, 3:6]
            else:
                # Unknown file variant
                msg = 'Unknown file format. Aborting..'
                raise Exception(msg)
            # end if

            return b, e
        # end func

        assert len(prefix), 'Error: invalid prefix. Aborting..'
        assert os.path.exists(outputPath), 'Error: invalid output path. Aborting..'

        if(len(list(self._filesDict.keys()))==0): return

        print('\nReading data files..')
        bxyz = None
        exyz = None
        for i, k in enumerate(self._filesDict.keys()):
            b, e = readData(self._filesDict[k])
            if(i == 0):
                if (b is not None): bxyz = b
                if (e is not None): exyz = e
            else:
                if(b is not None): bxyz = np.concatenate((bxyz, b), axis=0)
                if(e is not None): exyz = np.concatenate((exyz, e), axis=0)
            #end if
        #end for

        print('\nWriting output data files..')
        if(bxyz is not None):
            np.savetxt(os.path.join(outputPath, '%s.BX'%prefix), bxyz[:, 0], fmt='%.3f')
            np.savetxt(os.path.join(outputPath, '%s.BY' % prefix), bxyz[:, 1], fmt='%.3f')
            np.savetxt(os.path.join(outputPath, '%s.BZ' % prefix), bxyz[:, 2], fmt='%.3f')
        # end if

        if (exyz is not None):
            np.savetxt(os.path.join(outputPath, '%s.EX' % prefix), exyz[:, 0], fmt='%.3f')
            np.savetxt(os.path.join(outputPath, '%s.EY' % prefix), exyz[:, 1], fmt='%.3f')
            np.savetxt(os.path.join(outputPath, '%s.EZ' % prefix), exyz[:, 2], fmt='%.3f')
        # end if
    #end func
#end class

''' ========================================================
Setup Click interface
============================================================ '''
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('path', type=click.Path(exists=True))
@click.argument('output-path', type=click.Path(exists=True))
@click.option('--start-date-time', default='',
              help="Start date and time for selecting data files. Format should be YYYYMMDDHHMMSS; alternatively, "
                   "segments following a given trailing order term provided can be dropped. E.g. YYYY, YYYYMMDD, "
                   "etc. are valid.")
@click.option('--end-date-time', default='',
              help="End date and time for selecting data files, in the same format as --start-date-time.")
@click.option('--prefix', default='data', help="Prefix for output data file names.")
def process(path, output_path, start_date_time, end_date_time, prefix):
    '''
    PATH: Path to data files \n
    OUTPUT_PATH: Output folder

    Example: ./concatenate_input.py DATA0005/ /tmp/ --start-date-time 2014 --end-date-time 2017 --prefix='EX01'
    '''

    d = Data(path,
             startDateTime=start_date_time,
             endDateTime=end_date_time)
    d.ouput(prefix=prefix, outputPath=output_path)
# end func

if __name__ == "__main__":
    def test():
        d = Data('/home/rakib/work/ausLAMP/CT_workshop/testing2/DATA0005',
                 startDateTime='20170117',
                 endDateTime='20170118')
        d.ouput('a', '/tmp')
    # end func

    # Quick test
    if(0):
        test()
    else:
        process()
# end if
