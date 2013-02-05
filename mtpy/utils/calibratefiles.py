#!/usr/bin/env python
"""
This is a convenience script for the calibration of all files within a directory (non-recursive).

It needs the location of the directory and the location of the respective configuration file. If no output folder is specified, a subfolder 'calibrated' is set up within the input directory

"""

import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time


import mtpy.utils.exceptions as EX
reload(EX)

import mtpy.processing.calibration as C
reload(C)
import mtpy.processing.filehandling as FH
reload(FH)

angleaccuracy = 1.


def main():

    if len(sys.argv) < 3:
        raise EX.MTpyError_inputarguments('Need 2 arguments: <path to files> <config file>')


    pathname_raw = sys.argv[1] 
    directory = op.abspath(op.realpath(pathname_raw))

    configfilename_raw = sys.argv[2]
    configfile = op.abspath(op.realpath(configfilename_raw))


    if not op.isdir(directory):
        raise EX.MTpyError_inputarguments('Directory not existing: %s' % (directory))

    if not op.isfile(configfile):
        raise EX.MTpyError_inputarguments('Config file not existing: %s' % (configfile))


    try:
        outdir_raw = sys.argv[3]
        outdir = op.abspath(outdir_raw)
    except:
        outdir = op.join(directory,'calibrated')

    try:
        if not op.isdir(outdir):
            os.makedirs(outdir)
    except:
        raise EX.MTpyError_inputarguments('Output directory cannot be generated: %s' % (outdir))


    try:
        config_dir = FH.read_configfile(configfile)
    except:
        raise EX.MTpyError_config_file( 'Config file cannot be read: %s' % (configfile) )

    #----------------------------------------------------------------------------
    #to be improved later - rather rely on header lines than filenames!! :

    #select files by suffix, since header information is not necessarily present
    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz']

    oldwd = os.getcwd()
    os.chdir(directory)
    lo_allfiles = glob.glob('*.??')
    lo_allfiles = [op.abspath(i) for i in lo_allfiles]
    os.chdir(oldwd)

    lo_files = []

    for f in lo_allfiles:
        if f[-2:].lower() in components:
            lo_files.append(f)

    #check, if list of files is empty
    if len(lo_files) == 0:
        raise EX.MTpyError_inputarguments('Directory does not contain files to calibrate: %s' % (wd))
    #-------------------------------------------------------


    for f in lo_files:

        #find station 
        #try reading in a potentially existing header line
        try:
            F = open(f,'r')
            firstline = F.readline()
            F.close()
            firstlinesplit = firstline.strip().split()
            if firstlinesplit[0][0] == '#':
                #check for missing whitespace after commenting symbol #:
                if len(firstlinesplit[0]) == 1:
                    stationname = firstlinesplit[1].upper()
                    channel = firstlinesplit[2].lower()
                #otherwise take the rest of the first string as stationname    
                else:
                    stationname = firstlinesplit[0][1:].upper()
                    channel = firstlinesplit[1].lower()
            else:
                raise

        except:
            try:

                stationname = FH.EDL_get_stationname_fromfilename(f)
                channel = f[-2:].lower()
            except:
                print 'stationname or channel for file %s could not be determined - skipping file'%(f)
                continue
        

        #get configuration dictionary for this station
        try:
            stationdict = config_dir[stationname]
        except:
            print 'no entry for station %s found in configuration file %s skipping file'%(stationname, configfile )
            continue

        latitude = stationdict['latitude']
        longitude = stationdict['longitude']
        elevation = stationdict['elevation']

        field = channel[0]
        direction = channel[1]


        if field == 'e':
            #check North-South axis orientation
            if direction == 'x':
                angle = float(stationdict['e_xaxis_azimuth'])
                dipolelength = float(stationdict['e_xaxis_length'])
                if np.abs(180. - angle) < angleaccuracy:
                    #X-axis points southwards
                    dipolelength *= 1

                elif not np.abs(angle) < angleaccuracy:
                    print 'Configuration file error. E-field X-axis angle for station %s invalid: %f'%(stationname,angle)

            #check East-West axis orientation
            if direction == 'y':
                angle = float(stationdict['e_yaxis_azimuth'])
                dipolelength = float(stationdict['e_yaxis_length'])
                if np.abs(270. - angle) < angleaccuracy:
                    #X-axis points southwards
                    dipolelength *= 1
                elif not np.abs(90. - angle) < angleaccuracy:
                    print 'Configuration file error. E-field Y-axis angle for station %s invalid: %f'%(stationname,angle)

            logger = stationdict['e_logger_type']
            gain = stationdict['e_logger_gain']
            instrument = stationdict['e_instrument_type']
            instrument_amplification = stationdict['e_instrument_amplification']

        elif field == 'b':

            dipolelength = 1.
            logger = stationdict['b_logger_type']
            gain = stationdict['b_logger_gain']
            instrument = stationdict['b_instrument_type']
            instrument_amplification = stationdict['b_instrument_amplification']


        C.calibrate_file(f, outdir, instrument, logger, gain, dipolelength, stationname, channel, latitude, longitude, elevation,  offset = 0 )

if __name__=='__main__':
    main()