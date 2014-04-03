#!/usr/bin/env python

"""
/mtpy/utils/mseed.py

Functions for the conversion of raw time series: ASCII to/from miniSeed. 

The functionality is based on the "obspy.core" module from the ObsPy package.

Note: 
  Only single trace miniSeed are handled !


@UofA, 2013
(LK)

"""

#=================================================================

#import obspy.mseed as omseed
import numpy as np
from obspy.core import read, Trace, Stream, UTCDateTime
import os.path as op

#import pyrocko as pmseed

import mtpy.utils.exceptions as MTex
import mtpy.utils.filehandling as MTfh
import mtpy.utils.format as MTft
reload(MTft)
reload(MTex)
reload(MTfh)

#=================================================================



def convertfile_ts2miniseed(infile, outfile,channel=None, station = None, location=None, network = None):

    sta, cha, samplingrate, t_min, nsamples, unit, lat, lon, elev, data = MTfh.read_ts_file(infile)

    if station is None:
        station = sta.upper()
    if channel is None:
        channel = cha.upper()
    if location is None:
        location = ''
    if network is None:
        network = ''

    delta_t = 1./float(samplingrate)
    t0 = np.float64(t_min)

    outfilename = op.abspath(outfile)
    try:
        outfilename = writefile_obspy_singletrace(outfilename,station,channel,network,location, delta_t, t0, data)
    except:
        try:
            outfilename = writefile_pyrocko_singletrace(outfilename,station,channel,network,location, delta_t, t0, data)
        except:
            raise MTeX.MTpyError_inputarguments('ERROR - could not write minSeed file : {0}'.format(outfilename))

    return outfilename



def quadrupol_convertfile_miniseed2ts(infile, outfile,combine,invert):
    
    try:
        dummystream = read(infile)
    except:
        raise MTpyError_inputarguments('no valid miniSeed file')


    no_traces = len(dummystream)

    if no_traces < 4:
        print 'found only {0} traces in file - cannot combine all traces'


    lo_outfn = []
    lo_tuples = []
    for trace in range(no_traces):

        station, channel, location, network,  samplingrate, t0, nsamples, data =\
                                         readfile_obspy(infile,trace)

        channel = channel.lower()

        ts_tuple = [station.upper(), channel.lower(), samplingrate,t0, int(nsamples)]

        ts_tuple.append(data)

        lo_tuples.append(ts_tuple)

    if combine is True:

        lo_newtuples = []
        northtup = None
        southtup = None
        for tup in lo_tuples:
            if tup[1].lower()[-1] in ['n']:
                northtup = tup
                continue
            if tup[1].lower()[-1] in ['s']:
                southtup = tup
                continue
        if (northtup is not None) and (southtup is not None):
            newtup = [northtup[0], 'ex',northtup[2],northtup[3],northtup[4]]
            if invert is True:
                data = 0.5*(northtup[5] - southtup[5])
            else:
                data = 0.5*(northtup[5] + southtup[5])
            newtup.append(data)

            newtup = tuple(newtup)
        elif (northtup is not None):
            newtup = northtup
        else:
            newtup = southtup

        lo_newtuples.append(newtup)


        easttup = None
        westtup = None
        for tup in lo_tuples:
            if tup[1].lower()[-1] in ['e']:
                easttup = tup
                continue
            if tup[1].lower()[-1] in ['w']:
                westtup = tup
                continue

        if (easttup is not None) and (westtup is not None):
            newtup = [easttup[0], 'ey',easttup[2],easttup[3],easttup[4]]
            if invert is True:
                data = 0.5*(easttup[5] - westtup[5])
            else:
                data = 0.5*(easttup[5] + westtup[5])
            newtup.append(data)

            newtup = tuple(newtup)
        elif (easttup is not None):
            newtup = easttup
        else:
            newtup = westtup

        lo_newtuples.append(newtup)

    else:
        lo_newtuples = []
        for tup in lo_tuples:
            if tup[1].lower()[-1] in ['s','w']:
                if invert is True:
                    newtup = tuple([tup[:4],-tup[4]] )
                else:
                    newtup = tup
            else:
                newtup = tup
            lo_newtuples.append(newtup)

    for tup in lo_newtuples:

        outfilebase = op.splitext(op.abspath(outfile))[0]
        newoutfile = '{0}.{1}'.format(outfilebase,tup[1])


        outfilename = MTfh.write_ts_file_from_tuple(newoutfile,tup)
        #print 'wrote file {0}'.format(outfilename)
        
        lo_outfn.append(outfilename)

    return lo_outfn






def convertfile_miniseed2ts(infile, outfile, unit=None, lat = None, lon = None, elev = None):

    try:
        dummystream = read(infile)
    except: 
        raise MTpyError_inputarguments('infile is not miniSeed')

    no_traces = len(dummystream)
    lo_outfn = []


    for trace in range(no_traces):

        station, channel, location, network,  samplingrate, t0, nsamples, data =\
                                         readfile_obspy(infile,trace)

        channel = channel.lower()
        if channel[-1] == 'e':
            channel = channel[:-1]+ 'y'
        if channel[-1] == 'n':
            channel = channel[:-1] +'x'

        ts_tuple = [station.upper(), channel.lower(), samplingrate,t0, nsamples]


        if unit is not None:
            try:
                unit = unit.lower()
            except:
                unit = None

        # if unit is None:
        #     ts_tuple.append('unknown') 

        if lat is not None:
            try:
                lat = MTfT._assert_position_format('lat', lat)
            except:
                lat = None

        # if lat is None:
        #     ts_tuple.append(0.)


        if lon is not None:
            try:
                lon = MTfT._assert_position_format('lon', lon)
            except:
                lon = None

        # if lon is None:
        #     ts_tuple.append(0.)


        if elev is not None:
            try:
                elev = MTfT._assert_position_format('elev', elev)
            except:
                elev = None

        # if elev is None:
        #     ts_tuple.append(0.)


        ts_tuple.append(data)

        if outfile.lower().endswith('mseed'):
            outfilebase = op.splitext(op.abspath(outfile))[0]
            newoutfile = '{0}.{1}'.format(outfilebase,ts_tuple[1])
        else:
            newoutfile = outfile


        outfilename = MTfh.write_ts_file_from_tuple(newoutfile,tuple(ts_tuple))
        outfilename = newoutfile
        lo_outfn.append(outfilename)

    return lo_outfn



def readfile_obspy(infilename, trace = 0):
    
    infile = op.abspath(infilename)
    if not op.isfile(infile):
        raise MTeX.MTpyError_inputarguments('ERROR - miniSeed file not existing: {0}'.format(infile))

    try:
        ms_stream = read(infile)
    except:
        raise MTeX.MTpyError_inputarguments('ERROR - File is not a valid miniSed file: {0}'.format(infile))


    trace = ms_stream[trace]
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



def writefile_obspy_singletrace(outfilename,station,channel,network,location,
                                delta_t, t0, data):

    # Fill header attributes
    stats = {'network': network.upper(), 
             'station': station.upper(), 
             'location': location.upper(),
             'channel': channel.upper(), 
             'npts': len(data), 
             'sampling_rate': 1./delta_t, 
             'starttime' : t0}
    #define stream
    st = Stream([Trace(data=data, header=stats)])
    if not outfilename.lower().endswith('.mseed'):
        outfilename += '.mseed'
    
    

    #save to file
    outfilename = MTfh.make_unique_filename(outfilename)
    st.write(outfilename, 'MSEED')

    return outfilename