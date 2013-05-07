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
from obspy.core import read, Trace, Stream
import os.path as op

#import pyrocko as pmseed

import mtpy.utils.exceptions as EX
import mtpy.utils.filehandling as FH
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


def convertfile_miniseed2ts(infile, outfile, unit=None, lat = None, long = None, elev = None):

    
    outfilename = FH.write_ts_file_from_tuple(outfile,ts_tuple)

    return outfilename



def writefile_obspy_singletrace(outfilename,station,channel,network,location, delta_t, t0, data):

    # Fill header attributes
    stats = {'network': network.upper(), 'station': station.upper(), 'location': location.upper(),
         'channel': channel.upper(), 'npts': len(data), 'sampling_rate': 1./delta_t, 'starttime' : t0}
    #define stream
    st = Stream([Trace(data=data, header=stats)])
    if not outfilename.lower().endswith('mseed'):
        outfilename += '.mseed'
    #save to file
    i = 1
    while op.isfile(outfilename):
        outfilename = outfilename[:-6]+'_%i'%i+'.mseed'
        i += 1

    st.write(outfilename, format='MSEED')

    return outfilename