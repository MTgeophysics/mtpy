#!/usr/bin/env python
"""
This is a convenience script for the calibration of all (station-) files
within a directory (recursive/non-recursive).

Calibration includes the linear conversion from logger counts(mili volts)
into physical quantities, the conversion from voltage to electrtic field
strength for E-field components, and, with the appropriate option given,
the re-orientation of the data traces into true North(X)/East(Y) coordinates.

Required input: data files directory and the location of the respective
configuration file (survey meta-data file). If no output folder is specified,
 a subfolder 'calibrated' is set up within the input directory.

    2 mandatory arguments:
    - path to files
    - configuration file ('survey.cfg' style)

    4 optional arguments:
    - name of the output directory - cannot start with '-'
    - stationname - cannot start with '-'
    - flag '-R (or -r)', if the directory shall be searched for data recursively
    - flag '-O (or -o)', if the data shall be re-oriented to geographical N/E first

"""

import re
import sys
import os
import os.path as op


import mtpy.processing.calibration as MTcb
import mtpy.utils.filehandling as MTfh
import mtpy.utils.configfile as MTcf

#reload(MTfh)
#reload(MTcb)
#reload(MTcf)

# accuracy to which angles are determined
angleaccuracy = 1.


def main():

    if len(sys.argv) < 3:
        sys.exit('\nNeed at least 2 arguments:\n <path to files> \n '
                 '<config file> \n '
                 '[optional:<output dir>] \n [optional:<station>] \n '
                 '[optional:<recursive flag -R>] \n '
                 '[optional:<re-orientation flag -O]\n\n')

    outdir = None
    stationname = None
    recursive = False
    orientation = False

    if len(sys.argv) > 3:
        optionals = sys.argv[3:]
        for o in optionals:
            o = o.strip()
            if o[0] == '-':
                if o[1].lower() == 'r':
                    recursive = True
                elif o[1].lower() == 'o':
                    orientation = True
                continue
            elif outdir is None:
                outdir = o
                continue
            elif stationname is None:
                stationname = o.upper()
                continue
    pathname_raw = sys.argv[1]
    pathname = op.abspath(op.realpath(pathname_raw))

    if not op.isdir(pathname):
        sys.exit('Data file(s) path not existing: {0}'.format(pathname))

    configfilename_raw = sys.argv[2]
    configfile = op.abspath(
        op.realpath(
            op.join(
                os.curdir,
                configfilename_raw)))

    if not op.isfile(configfile):
        sys.exit('Config file not found: {0}'.format(configfile))

    if recursive is True:
        lo_dirs = [pathname]
        for i, j, k in os.walk(pathname):
            lof = [op.abspath(op.join(i, f)) for f in j]
            lo_dirs.extend(lof)
        pathname = list(set(lo_dirs))
    else:
        pathname = [pathname]

    #config_dict = MTcf.read_survey_configfile(configfile)
    try:
        config_dict = MTcf.read_survey_configfile(configfile)
        # done internally already
        # MTcf.validate_dict(config_dict)
    except:
        sys.exit(
            'Config file invalid or cannot be read: {0}'.format(configfile))

    #-------------------------------------------------------------------------

    # select files by header entries:
    components = ['ex', 'ey', 'bx', 'by', 'bz']
    lo_allfiles = []
    lo_allheaders = []
    lo_allstations = []
    for folder in pathname:
        wd = op.abspath(op.realpath(folder))
        if not op.isdir(wd):
            # print 'Directory not existing: %s' % (wd)
            lo_foldernames.remove(wd)
            continue
        dirfiles = [op.abspath(op.join(wd, i)) for i in os.listdir(wd)]
        for tmpfile in dirfiles:
            try:
                header = MTfh.read_ts_header(tmpfile)
                if header['channel'].lower() in components:
                    if stationname is not None:
                        if stationname.upper() != header['station'].upper():
                            continue
                    lo_allstations.append(header['station'].upper())
                    lo_allfiles.append(op.abspath(op.join(wd, tmpfile)))
                    lo_allheaders.append(header)
            except:
                continue

    lo_allstations = list(set(lo_allstations))

    # check, if list of files is empty
    if len(lo_allfiles) == 0:
        sys.exit('Directory(ies) do(es) not contain files to calibrate:'
                 ' {0}'.format(pathname))

    #-------------------------------------------------------
    # set up the directory structure for the output:

    # 1. generic calibration output directory
    cal_outdir = op.abspath(op.join(pathname[0], 'calibrated'))

    if outdir is not None:
        try:
            cal_outdir = op.abspath(op.join(os.curdir, outdir))
            if not op.isdir(cal_outdir):
                os.makedirs(cal_outdir)
                print('generated ', cal_outdir)
        except:
            print('Output directory cannot be generated: '\
                '{0} - using generic location'.format(cal_outdir))

            cal_outdir = op.abspath(op.join(pathname[0], 'calibrated'))
    try:
        if not op.isdir(cal_outdir):
            os.makedirs(cal_outdir)
    except:
        # this only comes up, if the generic location cannot be generated
        sys.exit(
            'Generic directory cannot be generated: {0}'.format(cal_outdir))

    print('\t Output directory ok: {0}\n'.format(cal_outdir))

    # if re-orientation is required, do it first:
    if orientation is True:
        print('\n\t....re-orient data first...\n')
        ori_outdir = op.abspath(op.join(cal_outdir, '../reoriented_tmp'))
        try:
            if not op.isdir(ori_outdir):
                os.makedirs(ori_outdir)
        except:
            # this only comes up, if the generic location cannot be generated
            sys.exit('Re-orientation directory cannot be generated:'
                     ' {0}'.format(ori_outdir))

        MTfh.reorient_files(lo_allfiles, configfile, lo_stations=lo_allstations,
                            outdir=ori_outdir)

        # change to calibration setup :
        outdir = cal_outdir
        new_inputdir = ori_outdir

        # file structure has changed, so the re-oriented files have to be read
        # again:
        components = ['ex', 'ey', 'bx', 'by', 'bz']
        lo_allfiles = []
        lo_allheaders = []
        lo_allstations = []
        dirfiles = [op.abspath(op.join(new_inputdir, i))
                    for i in os.listdir(new_inputdir)]
        for tmpfile in dirfiles:
            header = MTfh.read_ts_header(tmpfile)
            lo_allstations.append(header['station'].upper())
            lo_allfiles.append(tmpfile)
            lo_allheaders.append(header)

        lo_allstations = list(set(lo_allstations))

        # check, if list of files is empty
        if len(lo_allfiles) == 0:
            sys.exit('Directory(ies) do(es) not contain files to calibrate:'
                     ' {0}'.format(ori_outdir))

    #-------------------------------------------------
    # calibration

    lo_calibrated_files = []
    lo_calibrated_stations = []
    for file_idx, filename in enumerate(lo_allfiles):
        curr_station = lo_allheaders[file_idx]['station'].upper()
        if stationname is not None:
            if stationname.upper() != curr_station.upper():
                continue
        print('reading file {0}...'.format(filename))

        channel = lo_allheaders[file_idx]['channel']
        lo_calibrated_stations.append(curr_station)

        # get configuration dictionary for this station

        try:
            stationdict = config_dict[curr_station]
        except:
            print('no entry for station {0} found in configuration file'\
                ' {1} skipping file'.format(curr_station, configfile))
            continue

        latitude = float(stationdict['latitude'])
        longitude = float(stationdict['longitude'])
        elevation = float(stationdict['elevation'])

        field = channel[0]
        direction = channel[1]

        station_type = stationdict['station_type']

        if field == 'e':
            if station_type == 'b':
                continue
            # check North-South axis orientation
            if direction == 'x':
                #angle = float(stationdict['e_xaxis_azimuth'])
                dipolelength = float(stationdict['e_xaxis_length'])

            # check East-West axis orientation
            if direction == 'y':
                #angle = float(stationdict['e_yaxis_azimuth'])
                dipolelength = float(stationdict['e_yaxis_length'])

            logger = stationdict['e_logger_type']
            gain = float(stationdict['e_logger_gain'])
            instrument = stationdict.get('e_instrument_type', 'electrodes')
            instrument_amplification = float(
                stationdict['e_instrument_amplification'])

        elif field == 'b':
            if station_type == 'e':
                continue

            dipolelength = 1.
            logger = stationdict['b_logger_type']
            gain = float(stationdict['b_logger_gain'])
            instrument = stationdict.get('b_instrument_type', 'coils')
            instrument_amplification = float(
                stationdict['b_instrument_amplification'])

        MTcb.calibrate_file(filename, outdir, instrument, instrument_amplification,
                            logger, gain, dipolelength, curr_station, channel,
                            latitude, longitude, elevation, offset=0)
        # print 'calibrated file {0},{1}'.format(outdir, filename)
        lo_calibrated_files.append(filename)

    lo_calibrated_stations = list(set(lo_calibrated_stations))
    if len(lo_calibrated_files) == 0:
        if stationname is not None:
            print('No files found for station {0}'.format(stationname))
            return
        else:
            print('No files found for stations {0}'.format(lo_allstations))

    print('{0} files calibrated for stations'\
        ' {1}'.format(len(lo_calibrated_files), lo_calibrated_stations))

if __name__ == '__main__':
    main()
