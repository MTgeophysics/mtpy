#!/usr/bin/env python

"""


"""


#import gc
from numpy import *
import sys
import os
import os.path as op

from pyrocko import trace, util, io

#import pdb

channel_dict = {'n': 0, 'e': 1, 's': 2, 'w': 3}

longchannels = ['north', 'east', 'south', 'west']


def main():

    if len(sys.argv) < 5:
        sys.exit('\n\tERROR - need 4 arguments as input: \n <input dir> <output dir>'
                 '<stationname> <prefix>\n')
    print()

    indir = sys.argv[1]
    indir = op.abspath(op.join(os.curdir, indir))

    outdir = sys.argv[2]
    outdir = op.abspath(op.join(os.curdir, outdir))

    stationname = sys.argv[3].upper()

    fullprefix = sys.argv[4]
    try:
        prefix = fullprefix.split('.')[0]
    except:
        prefix = fullprefix

    try:
        t0, samplingrate = read_timestampfile(indir, prefix)
    except:
        sys.exit(
            '\n ERROR - cannot find timestamp file for given prefix: {0}\n'.format(prefix))

    if not op.isdir(outdir):
        os.makedirs(outdir)

    lo_infiles = os.listdir(indir)
    lo_infiles = [
        i for i in lo_infiles if i.lower().startswith(
            fullprefix.lower())]

    for fn in lo_infiles:
        ext = op.splitext(fn)[-1][1:].lower()
        if 1:
            if ext[0] in 'nesw':
                chan = ext[0]
                infile = op.join(indir, fn)
                outfn = makeoutfilename(outdir, stationname, prefix, chan)
                print(outfn)
                run(infile, outfn, stationname, t0, samplingrate, chan=chan)
                print()
        # except:
        #     continue


def makeoutfilename(outdir, stationname, prefix, chan):

    longchan = longchannels[channel_dict[chan]]

    outfilename = '%s_%s.%s.mseed' % (stationname.upper(), prefix, longchan)
    outfn = op.join(outdir, outfilename)

    if op.exists(outfn):
        file_exists = True
        number = 1
        while file_exists:
            outfilename = '%s_%s.%s.%i.mseed' % (
                stationname.upper(), prefix, longchan, number)
            outfn = op.join(outdir, outfilename)
            if not op.exists(outfn):
                file_exists = False
            number += 1

    return outfn


def read_timestampfile(indir, prefix):

    fn = prefix + '.timestamps'
    fn = op.join(indir, fn)
    t0 = None
    sampling = None
    in_dict = {}
    Fin = open(fn)
    for line in Fin:
        line = line.strip()
        try:
            in_dict[line.split(':')[0].strip()] = line.split(':')[1].strip()
        except:
            continue

    for k, v in list(in_dict.items()):
        if k.lower().startswith('sampl'):
            sampling = float(v)
            if sampling % 1 == 0:
                sampling = int(sampling)

        if k.lower().startswith('first'):
            t0 = float64(v)

    return t0, sampling


def run(infile, outfn, station, t0, samplingrate, chan='n', nw='', loc=''):

    channels = 'NESW'

    # read in data:
    data = []
    print('\t reading file {0} ... '.format(infile))
    Fin = open(infile)
    for line in Fin:
        if line.strip()[0] == '#':
            continue
        data.append([int(float(i)) for i in line.strip().split()])
    Fin.close()
    data = array(data, dtype=int32)
    # print '...done!'

    stationname = station
    deltat = 1. / samplingrate

    lo_traces = []

    lo_chans = []
    try:
        lo_chans = list(chan.lower())
    except:
        lo_chans = ['n']

    print('\t building MiniSeed trace object(s)...')
    if len(lo_chans) > 1:
        try:
            for c in lo_chans:
                location = c.upper()
                idx_ch = channel_dict[c]
                if idx_ch in [0, 2]:
                    channel = 'NS'
                else:
                    channel = 'EW'

            print('station {1} - channel {0} - location {2}'.format(channel, stationname, location))
            if idx_ch in [0, 1]:
                lo_traces.append(trace.Trace(station=stationname, channel=channel,
                                             location=location, deltat=deltat, tmin=t0, ydata=data[:, idx_ch]))
            else:
                # correct for polarity of 'redundant' channels
                lo_traces.append(trace.Trace(station=stationname, channel=channel,
                                             location=location, deltat=deltat, tmin=t0, ydata=-data[:, idx_ch]))
        except:
            lo_chans = [chan.lower()[0]]
            lo_traces = []

    if len(lo_chans) == 1:
        idx_ch = channel_dict[lo_chans[0]]
        location = channels[idx_ch]
        if idx_ch in [0, 2]:
            channel = 'NS'
        else:
            channel = 'EW'

        if idx_ch in [2, 3]:
            data *= -1

        lo_traces = [trace.Trace(station=stationname, channel=channel,
                                 location=location, deltat=deltat, tmin=t0, ydata=data)]

    # pdb.set_trace()
    # print '...done!'

    print('\t writing file %s ... ' % outfn)
    io.save(lo_traces, outfn)
    print('\t ...done')
    return


if __name__ == '__main__':
    main()
