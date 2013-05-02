#!/usr/bin/env python
"""
This is a convenience script for the calibration of all (station-) files within a directory (recursive/non-recursive).

Calibration includes the linear conversion from logger counts(mili volts) into physical quantities, the conversion from voltage to electrtic field strength for Efield components, and the re-orientation of the data traces into true North(X)/East(Y) coordinates.

Required input: data files directory and the location of the respective configuration file (survey meta-data file). If no output folder is specified, a subfolder 'calibrated' is set up within the input directory.

    2 mandatory arguments: 
    - path to files 
    - configuration file ('survey.cfg' style)

    3 optional arguments:
    - name of the output directory - cannot start with '-' 
    - stationname - cannot start with '-' 
    - flag '-R (or -r)', if the directory shall be searched for data recursively 

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
import mtpy.processing.calibration as CAL
import mtpy.processing.general as GEN
import mtpy.utils.filehandling as FH

reload(FH)
reload(EX)
reload(CAL)
reload(GEN)

angleaccuracy = 1.


def main():

    if len(sys.argv) < 3:
        raise EX.MTpyError_inputarguments('Need at least 2 arguments: <path to files> <config file>  [<output dir>] [<stationname>] [<recursive flag -R>]')
    outdir = None
    stationname = None
    recursive = False

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

    pathname_raw = sys.argv[1] 
    pathname = op.abspath(op.realpath(pathname_raw))

    if not op.isdir(pathname):
        raise EX.MTpyError_inputarguments('Data file(s) path not existing: {0}'.format(pathname))

    configfilename_raw = sys.argv[2]
    configfile = op.abspath(op.realpath(op.join(os.curdir,configfilename_raw)))

    if not op.isfile(configfile):
        raise EX.MTpyError_inputarguments('Config file not found: {0}'.format(configfile))


    if recursive is True:
        lo_files = []
        for i,j,k in os.walk(pathname):
            lof = [op.abspath(op.join(i,f)) for f in j]
            lo_files.extend(lof)
        pathname = list(set(lo_files))
    else:
        pathname = [pathname]

    if outdir is not None:
        try:
            outdir = op.abspath(op.join(os.curdir,outdir))
            if not os.isdir(outdir):
                os.makedirs(outdir)
        except:
            outdir = op.join(pathname[0],'calibrated')
            if not os.isdir(outdir):
                try:
                    os.makedirs(outdir)
                except:
                    EX.MTpyError_inputarguments('Output directory cannot be generated: %s' % (outdir))

    try:
        config_dict = FH.read_configfile(configfile)
    except:
        raise EX.MTpyError_config_file( 'Config file cannot be read: %s' % (configfile) )

    #----------------------------------------------------------------------------

    #select files by heder entries:
    components = ['ex', 'ey', 'bx', 'by', 'bz']
    lo_allfiles = []
    lo_allheaders = []
    for folder in pathname:
        wd = op.abspath(op.realpath(folder)) 
        if not op.isdir(wd):
            #print 'Directory not existing: %s' % (wd)
            lo_foldernames.remove(wd)
            continue    
        dirfiles = os.listdir(wd)
        for tmpfile in dirfiles:
            header = FH.read_ts_header(tmpfile)
            if header['channel'].lower() in components:
                lo_allfiles.append(op.abspath(op.join(wd,tmpfile)))
                lo_allheaders.append(header)


    #check, if list of files is empty
    if len(lo_allfiles) == 0:
        raise EX.MTpyError_inputarguments('Directory(ies) do(es) not contain files to calibrate: {0}'.format(pathname))
    #-------------------------------------------------------



        # #find station 
        # #try reading in a potentially existing header line
        # try:
        #     F = open(filename,'r')
        #     firstline = F.readline()
        #     F.close()
        #     firstlinesplit = firstline.strip().split()
        #     if firstlinesplit[0][0] == '#':
        #         #check for missing whitespace after commenting symbol #:
        #         if len(firstlinesplit[0]) == 1:
        #             stationname = firstlinesplit[1].upper()
        #             channel = firstlinesplit[2].lower()
        #         #otherwise take the rest of the first string as stationname    
        #         else:
        #             stationname = firstlinesplit[0][1:].upper()
        #             channel = firstlinesplit[1].lower()
        #     else:
        #         raise

        # except:
        #     try:

        #         stationname = FH.EDL_get_stationname_fromfilename(f)
        #         channel = filename[-2:].lower()
        #     except:
        #         print 'stationname or channel for file {0} could not be determined - skipping file'.format(filename)
        #         continue

    for file_idx, filename in enumerate(lo_allfiles):
        
        stationname = lo_allheaders[file_idx]['station']
        channel = lo_allheaders[file_idx]['channel']

        #get configuration dictionary for this station
        try:
            stationdict = config_dir[stationname]
        except:
            print 'no entry for station {0} found in configuration file {1} skipping file'.format(stationname, configfile )
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

            #check East-West axis orientation
            if direction == 'y':
                angle = float(stationdict['e_yaxis_azimuth'])
                dipolelength = float(stationdict['e_yaxis_length'])

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


        CAL.calibrate_file(filename, outdir, instrument, logger, gain, dipolelength, stationname, channel, latitude, longitude, elevation,  offset = 0 )

TODO : re-orientate files !!

        GEN.correct4sensor_orientation(x_values, y_values, x_sensor_angle = 0 , y_sensor_angle = 90):
if __name__=='__main__':
    main()