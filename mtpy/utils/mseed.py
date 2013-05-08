#!/usr/bin/env python

"""
/mtpy/utils/mseed.py

This modules contains functions for the conversion of raw time series: ASCII to/from miniSeed. 

The functionality is based on obspy.mseed/pyrocko


@UofA, 2013
(LK)

"""

#=================================================================

#import obspy.mseed as omseed
from obspy.core import read, Trace, Stream, UTCDateTime
import os.path as op

#import pyrocko as pmseed

import mtpy.utils.exceptions as EX
import mtpy.utils.filehandling as FH
import mtpy.utils.format as FT
reload(FT)
reload(EX)
reload(FH)

#=================================================================



def convertfile_ts2miniseed(infile, outfile,channel=None, station = None, location=None, network = None):

    sta, cha, samplingrate, t_min, nsamples, unit, lat, lon, elev, data = FH.read_ts_file(infile)

    if station is None:
        station = sta.upper()
    if channel is None:
        channel = cha.upper()
    if location is None:
        location = ''
    if network is None:
        network = ''

    delta_t = 1./float(samplingrate)
    t0 = float(t_min)

    outfilename = op.abspath(outfile)
    try:
        outfilename = writefile_obspy_singletrace(outfilename,station,channel,network,location, delta_t, t0, data)
    except:
        try:
            outfilename = writefile_pyrocko_singletrace(outfilename,station,channel,network,location, delta_t, t0, data)
        except:
            raise EX.MTpyError_inputarguments('ERROR - could not write minSeed file : {0}'.format(outfilename))

    return outfilename


def convertfile_miniseed2ts(infile, outfile, unit=None, lat = None, lon = None, elev = None):

    station, channel, location, network,  samplingrate, t0, nsamples, data = readfile_obspy_singletrace(infile)

    ts_tuple = [station.upper(), channel.lower(), samplingrate,t0, nsamples]

    if unit is not None:
        try:
            unit = unit.lower()
        except:
            unit = None

    if unit is None:
        ts_tuple.append('unknown') 

    if lat is not None:
        try:
            lat = FT._assert_position_format('lat', lat)
        except:
            lat = None

    if lat is None:
        ts_tuple.append(0.)


    if lon is not None:
        try:
            lon = FT._assert_position_format('lon', lon)
        except:
            lon = None

    if lon is None:
        ts_tuple.append(0.)


    if elev is not None:
        try:
            elev = FT._assert_position_format('elev', elev)
        except:
            elev = None

    if elev is None:
        ts_tuple.append(0.)

    ts_tuple.append(data)

    print data[:10]

    outfilename = FH.write_ts_file_from_tuple(outfile,tuple(ts_tuple))

    return outfilename



def readfile_obspy_singletrace(infilename):
    
    infile = op.abspath(infilename)
    if not op.isfile(infile):
        raise EX.MTpyError_inputarguments('ERROR - miniSeed file not existing: {0}'.format(infile))

    try:
        ms_stream = read(infile)
    except:
        EX.MTpyError_inputarguments('ERROR - File is not a valid miniSed file: {0}'.format(infile))

    trace = ms_stream[0]
    stats = trace.stats

    t0 = stats['starttime'].timestamp
    station = stats['station']
    channel = stats['channel']
    location = stats['location']
    network = stats['network']
    samplingrate = 1./float(stats['delta'])
    nsamples = trace.count()

    data = trace.data

    return station, channel, location, network,  samplingrate, t0, nsamples, data



def writefile_obspy_singletrace(outfilename,station,channel,network,location, delta_t, t0, data):

    # Fill header attributes
    stats = {'network': network.upper(), 'station': station.upper(), 'location': location.upper(),
         'channel': channel.upper(), 'npts': len(data), 'sampling_rate': 1./delta_t, 'starttime' : t0}
    #define stream
    st = Stream([Trace(data=data, header=stats)])
    if not outfilename.lower().endswith('.mseed'):
        outfilename += '.mseed'
    

    #save to file
    outfilename = FH.make_unique_filename(outfilename)
    st.write(outfilename, format='MSEED')

    return outfilename