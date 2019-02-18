#!/usr/bin/env python

import numpy as np
import re
import os
import sys
import os.path as op
import shutil
import calendar
import time
import subprocess

#import ipdb

channel_dict = {'n': 0, 'e': 1, 's': 2, 'w': 3}

scriptbase = '/data/software/SymRes_scripts/May2014'

scriptbase = op.abspath(op.join(os.curdir, scriptbase))

# output data block length (in hours)
outputlength = 24

samplingrate = 1

samplingrate_keys = {'1': 0, '500': 8, '200': 5}


cleanup = True
#cleanup = False


def main():
    """

    """

    if len(sys.argv) < 8:
        print("""\n use 7 arguments:
        <data directory sta A> <data directory sta B> <data directory sta C>
        <channelsA> <channelsB> <channelsC>
        <start time yymmdd-HHhMM> [<destination folder>]

        """)
        sys.exit()

    # deal with command line arguments
    dataDirA = sys.argv[1]
    dataDirB = sys.argv[2]
    dataDirC = sys.argv[3]
    channelsAin = sys.argv[4]
    channelsBin = sys.argv[5]
    channelsCin = sys.argv[6]
    timestring = sys.argv[7]

    # check, if output directory has been given
    if len(sys.argv) > 8:
        dst_raw = sys.argv[8]
    try:
        dst = op.join(op.abspath(os.curdir), dst_raw)
    except:
        dst = op.join(
            op.abspath(
                os.curdir),
            'birrpdata_%dh' %
            (int(outputlength)))

    print()
    if not op.isdir(dst):
        os.makedirs(dst)

    if (not op.isdir(dataDirA)):
        sys.exit("\nERROR - input data A directory not found \n")
    if (not op.isdir(dataDirB)):
        sys.exit("\nERROR - input data B directory not found \n")
    if (not op.isdir(dataDirB)):
        sys.exit("\nERROR - input data C directory not found \n")

    allchannels = 'nesw'
    channelsA = []
    tmp_chs = channelsAin.split(',')
    try:
        if len(tmp_chs) < 2:
            raise
        for i in tmp_chs:
            ch = i.lower()
            if not ch in allchannels:
                raise
            else:
                channelsA.append(channel_dict[ch])
        if len(channelsA) == 0:
            raise
    except:
        channelsA = [0, 1, 2, 3]

    channelsB = []
    tmp_chs = channelsBin.split(',')
    try:
        if len(tmp_chs) < 2:
            raise
        for i in tmp_chs:
            ch = i.lower()
            if not ch in allchannels:
                raise
            else:
                channelsB.append(channel_dict[ch])
        if len(channelsB) == 0:
            raise
    except:
        channelsB = [0, 1]

    channelsC = []
    tmp_chs = channelsCin.split(',')
    try:
        if len(tmp_chs) < 2:
            raise
        for i in tmp_chs:
            ch = i.lower()
            if not ch in allchannels:
                raise
            else:
                channelsC.append(channel_dict[ch])
        if len(channelsC) == 0:
            raise
    except:
        channelsC = [0, 1]

    try:
        timestamp = time.strptime(timestring, "%y%m%d-%Hh%M")
    except:
        sys.exit(
            "\n ERROR - unknown time - check format: yymmdd-HHhMM (e.g. 140123-06h00)\n")

    # call the actuall processing:
    run(dataDirA, dataDirB, dataDirC, channelsA, channelsB, channelsC, timestamp,
        outputlength, samplingrate, dst, scriptbase)


def run(dataDirA, dataDirB, dataDirC, channelsA, channelsB, channelsC, timestamp, duration,
        samplingrate, dst, SymResScriptbase, keepfiles=True):

    sampling_mus = int(1. / samplingrate * 1e6)

    starttime = calendar.timegm(timestamp)
    endtime = starttime + (duration * 3600)

    # print starttime,endtime

    lo_subdirsA, firstsampleA, endA = findDataSubDirs(
        dataDirA, starttime, endtime)
    lo_subdirsB, firstsampleB, endB = findDataSubDirs(
        dataDirB, starttime, endtime)
    lo_subdirsC, firstsampleC, endC = findDataSubDirs(
        dataDirC, starttime, endtime)

    firstsample = max([firstsampleA, firstsampleB, firstsampleC, starttime])
    finalsample = min([endA, endB, endC, endtime]) - 1e-6 * sampling_mus

    # print firstsampleA,firstsampleB,firstsampleC,starttime
    # print firstsample

    # print endA,endB,endC,endtime
    # print finalsample

    lo_lo_subdirs = [lo_subdirsA, lo_subdirsB, lo_subdirsC]

    lo_dirs = [dataDirA, dataDirB, dataDirC]
    # update file lists using updated start and end times

    basedir = op.abspath(os.curdir)

    stationlist = ['A', 'B', 'C']

    # print lo_lo_subdirs
    # sys.exit()
    print()
    # loop over 3 stations
    for idx_station, station in enumerate(stationlist):

        os.chdir(lo_dirs[idx_station])
        print('\t\tstation', station, op.abspath(os.curdir))
        print()

        lo_subdirs = lo_lo_subdirs[idx_station]

        current_last_sample = firstsample

        # all the following actions require unpacking and SymRes-processing
        # of the respective raw data files.
        # The option 'keepfiles=True' allows to keep the unpacked and interpolated
        # ASCII data file. This allows for faster re-processing if needed.

        # get data from the first directory - from sample matching the
        # starttime
        outdata, current_last_sample = evaluate_first_directory(lo_subdirs[0], samplingrate,
                                                                firstsample, SymResScriptbase)

        # print outdata.shape,current_last_sample

        # add data from following directories - including interpolation within
        # gaps
        for i in lo_subdirs[1:-1]:
            tmp_data, current_last_sample = evaluate_intermediate_directory(i, samplingrate,
                                                                            outdata[-1], current_last_sample, SymResScriptbase)
            outdata = np.append(outdata, tmp_data, axis=0)
            # print outdata.shape,current_last_sample

        # add data from the last directory - incl. interpolation - until
        # endtime stamp
        tmp_data = evaluate_last_directory(lo_subdirs[-1], samplingrate, outdata[-1],
                                           current_last_sample, finalsample, SymResScriptbase)
        outdata = np.append(outdata, tmp_data, axis=0)
        # print outdata.shape

        # write data to file
        outfn = 'data%s' % (station)
        outfn = op.join(dst, outfn)
        print('\n    Writing data to file: {0} ...'.format(outfn))

        Fout = open(outfn, 'w')

        idx = 0
        for line in outdata:
            idx += 1
            try:
                for item in line:
                    Fout.write('%d    ' % int(item))
            except:
                print(idx, line)
                sys.exit()
            Fout.write('\n')

        Fout.close()
        print('\t ... Done !')

        os.chdir(basedir)
        print()
        print('\t' + 72 * '-')
        print()

    # sys.exit()
    # split up files within the destination folder
    prepare_data(
        dst,
        channelsA,
        channelsB,
        channelsC,
        firstsample,
        finalsample)

    if cleanup is False:
        print('\n No cleanup!\n\n')
        sys.exit()

    # loop over all dis and subdirs for cleaning up:
    print()
    print('\n\t\tcleaning up files ...\n')
    for idx_station, station in enumerate(stationlist):
        print(station)
        os.chdir(lo_dirs[idx_station])
        stationdir = op.abspath('.')
        lo_subdirs = lo_lo_subdirs[idx_station]
        for subdir in lo_subdirs:
            os.chdir(subdir)
            print('\t', op.abspath('.'))
            cleanupfiles(keepfiles)
            os.chdir(stationdir)
            # print op.abspath('.')
        print()
        os.chdir(basedir)

    # print '\n\t removing temporary data files ...\n'
    # remove_datafiles(dst)

    print('\n\n\t\t!! DONE !!\n')
    print(80 * '=' + '\n')


def remove_datafiles(dst):

    cwd = op.abspath(os.curdir)
    os.chdir(dst)
    lo_files = os.listdir('.')
    lo_files = [i for i in lo_files if i.startswith('data')]
    for f in lo_files:
        os.remove(f)
    os.chdir(cwd)


def directoryname2timestamp(dirname):

    return calendar.timegm(time.strptime(dirname, "%Y-%m-%d-at-%H-%M-%S.bz"))


def findDataSubDirs(dataDir, starttime, endtime):

    lo_subdirs = sorted(os.listdir(dataDir))
    lo_subdirs = [i for i in lo_subdirs if i.endswith('.bz')]

    lo_starttimes = []
    lo_dirs = []
    for idx, curdir in enumerate(lo_subdirs):
        dir_start = directoryname2timestamp(curdir)
        if starttime <= dir_start < endtime:
            # if it's the first of interest:
            if len(lo_dirs) == 0:
                # check, if the directory before existed:
                if idx > 0:
                    # check, if the one before ended after the starttime:
                    tempdir = lo_subdirs[idx - 1]
                    tempstart = directoryname2timestamp(tempdir)
                    tempend = find_directory_endtime(
                        op.join(dataDir, tempdir), tempstart)
                    if tempend > starttime:
                        lo_dirs.append(tempdir)
            lo_dirs.append(curdir)

    # now find the overall endtime in the same way:
    firststart = directoryname2timestamp(lo_dirs[0])
    laststart = directoryname2timestamp(lo_dirs[-1])
    overall_endtime = find_directory_endtime(
        op.join(dataDir, lo_dirs[-1]), laststart)

    # directory names are only exact to the minute, so make sure that no
    # data are request that are not in here:
    firststart += 60

    return lo_dirs, firststart, overall_endtime


def find_directory_endtime(directory, starttime):

    lo_packedpaks = os.listdir(directory)
    lo_packedpaks = [i for i in lo_packedpaks if i.endswith('.pak.bz2')]
    lo_paks = os.listdir(directory)
    lo_paks = [i for i in lo_paks if i.endswith('.pak')]
    # found number of 10min blocks....now find the end time of
    # this subdirectory:
    num_10minblocks = max(len(lo_paks), len(lo_packedpaks))

    tempend = starttime + 600 * num_10minblocks

    return tempend


def unpackandinterpolate(SymResScriptbase, sammplingrate):

    sampling_key = samplingrate_keys['%d' % (int(samplingrate))]

    print('unpacking files in directory {0}'.format(op.abspath(os.curdir)))

    lo_files = sorted(os.listdir('.'))
    for f in lo_files:
        if f.lower().endswith('bz2'):
            try:
                os.system("bunzip2 -k -f {0}".format(f))
            except:
                continue
            print('unpacked {0}'.format(f))

    print('converting to PakBin file')

    print(SymResScriptbase)
    pak2bin_string = op.join(SymResScriptbase, "Pak2Bin")
    os.system('"{0}" 00000000.pak'.format(pak2bin_string))

    print('interpolating to {0} Hz'.format(samplingrate))
    interp_string = op.join(SymResScriptbase, "InterpAdl")
    os.system(
        '"{0}" Pak2Bin-300-Data.bin o{1}'.format(interp_string, sampling_key))


def evaluate_first_directory(
        directory, samplingrate, firstsample, SymResScriptbase):

    ascii_filename = 'QELPakData.interp%d.adl' % (int(samplingrate))
    fullname = op.join(directory, ascii_filename)

    cwd = op.abspath(os.curdir)
    os.chdir(directory)
    if not op.isfile(ascii_filename):
        unpackandinterpolate(SymResScriptbase, samplingrate)

    data = []

    Fin = open(ascii_filename)
    idx = 0
    for line in Fin:
        t = np.float64(line.strip().split()[-1])
        if t < firstsample:
            continue
        if idx == 0:
            t0 = t
            idx += 1

        data.append([int(float(i)) for i in line.strip().split()[:-1]])

    # ipdb.set_trace()

    stepsize = 1. / samplingrate
    lastsample = t0 + (len(data) - 1) * stepsize
    # print samplingrate,stepsize,t0,firstsample,len(data),lastsample

    print('\t\tfirst directory: {3} --\n\t{0} data points from {1:.3f} up to t= {2:.3f}\n'.format(
        len(np.array(data)), t0, lastsample, directory))

    os.chdir(cwd)

    return np.array(data), lastsample


def cleanupfiles(keepfiles):

    allfiles = sorted(os.listdir('.'))

    # wipe unpacked bz files, only if original is still present
    bzfiles = [i for i in allfiles if i.endswith('.bz2')]
    if len(bzfiles) > 0:
        files2kill2 = [i for i in allfiles if i.endswith('.pak')]
        for f in files2kill2:
            os.remove(f)
    # delete unpack-/interpolation auxiliary files
    files2kill = [i for i in allfiles if i.startswith('Pak2')]
    for f in files2kill:
        os.remove(f)

    # remove interpolated ASCII data, if flagged to do so
    if keepfiles is False:
        files2kill3 = [i for i in allfiles if i.startswith('QEL')]
        for f in files2kill3:
            os.remove(f)


def evaluate_intermediate_directory(
        directory, samplingrate, olddata, oldend, SymResScriptbase):

    ascii_filename = 'QELPakData.interp%d.adl' % (int(samplingrate))
    fullname = op.join(directory, ascii_filename)

    cwd = op.abspath(os.curdir)
    os.chdir(directory)

    # print op.abspath(os.curdir)#directory,oldend,samplingrate,keepfiles

    if not op.isfile(ascii_filename):
        unpackandinterpolate(SymResScriptbase, samplingrate)

    data = []

    linenumber = 0
    Fin = open(ascii_filename)
    firstline = Fin.readline().strip().split()
    linenumber += 1
    firstvalues = [int(float(i)) for i in firstline[:-1]]

    t0 = np.float64(firstline[-1])

    gap = t0 - oldend
    no_gapsamples = int(round(gap * samplingrate) - 1)
    print(30 * ' ' + '--->')
    print('\tneed {0} samples to bridge gap between {1:.3f} and {2:.3f}'.format(no_gapsamples,
                                                                                oldend, t0))

    # linspace includes start- and endpoint
    ch1 = np.linspace(olddata[0], firstvalues[0], no_gapsamples + 2)[1:]
    ch2 = np.linspace(olddata[1], firstvalues[1], no_gapsamples + 2)[1:]
    ch3 = np.linspace(olddata[2], firstvalues[2], no_gapsamples + 2)[1:]
    ch4 = np.linspace(olddata[3], firstvalues[3], no_gapsamples + 2)[1:]

    for i in range(no_gapsamples + 1):
        data.append([int(round(ch1[i], 0)), int(round(ch2[i])),
                     int(round(ch3[i])), int(round(ch4[i]))])

    for line in Fin.readlines():
        data.append([int(float(i)) for i in line.strip().split()[:-1]])
        linenumber += 1

    Fin.close()

    stepsize = 1. / samplingrate
    lastsample = t0 + (linenumber - 1) * stepsize

    print('\t\tintermediate directory: {2} -- \n\t{0} data points in file up to t= {1:.3f}\n'.format(linenumber, lastsample, directory))

    os.chdir(cwd)
    return np.array(data), lastsample


def evaluate_last_directory(directory, samplingrate,
                            olddata, oldend, finalsample, SymResScriptbase):

    ascii_filename = 'QELPakData.interp%d.adl' % (int(samplingrate))
    fullname = op.join(directory, ascii_filename)

    cwd = op.abspath(os.curdir)
    os.chdir(directory)

    # print op.abspath(os.curdir)#directory,oldend,samplingrate,keepfiles

    if not op.isfile(ascii_filename):
        unpackandinterpolate(SymResScriptbase, samplingrate)

    data = []
    linenumber = 0
    Fin = open(ascii_filename)
    firstline = Fin.readline().strip().split()
    linenumber += 1
    firstvalues = [int(float(i)) for i in firstline[:-1]]

    t0 = np.float64(firstline[-1])

    gap = t0 - oldend
    no_gapsamples = int(round(gap * samplingrate) - 1)
    print(30 * ' ' + '--->')
    print('\tneed {0} samples to bridge gap between {1:.3f} and {2:.3f}'.format(
        no_gapsamples, oldend, t0))
    print('\t\tlast directory: {0} -- '.format(directory))

    # linspace includes start- and endpoint
    ch1 = np.linspace(olddata[0], firstvalues[0], no_gapsamples + 2)[1:]
    ch2 = np.linspace(olddata[1], firstvalues[1], no_gapsamples + 2)[1:]
    ch3 = np.linspace(olddata[2], firstvalues[2], no_gapsamples + 2)[1:]
    ch4 = np.linspace(olddata[3], firstvalues[3], no_gapsamples + 2)[1:]

    for i in range(no_gapsamples + 1):
        data.append([int(round(ch1[i], 0)), int(round(ch2[i])),
                     int(round(ch3[i])), int(round(ch4[i]))])

    for line in Fin.readlines():
        t = np.float64(line.strip().split()[-1])
        if t > finalsample:
            continue
        data.append([int(float(i)) for i in line.strip().split()[:-1]])
        linenumber += 1

    print('\t{0} data points in file up to t= {1:.3f}'.format(linenumber, finalsample))

    Fin.close()

    return np.array(data)


def prepare_data(dst, channelsA, channelsB,
                 channelsC, firstsample, lastsample):

    fileA = 'dataA'
    fileB = 'dataB'
    fileC = 'dataC'

    channels = ['north', 'east', 'south', 'west']
    lo_stations = ['A', 'B', 'C']
    lo_lo_channels = [channelsA, channelsB, channelsC]

    print('\t' + 72 * '-' + '\n')
    print('\t\twriting data to single column files...')

    basedir = op.abspath(os.curdir)

    os.chdir(dst)

    tA = time.gmtime(firstsample)
    filebase = '{0:02d}{1:02d}{2:02d}-{3:02d}h{4:02d}'.format(
        tA[0] % 100, tA[1], tA[2], tA[3], tA[4])

    for idx_sta, station in enumerate(lo_stations):
        infile = 'data' + station
        for idx_c, ch in enumerate(channels):
            stationchannels = lo_lo_channels[idx_sta]
            if idx_c in stationchannels:
                # print 'pid ({0})'.format(os.getpid()), idx_c,ch

                fn_out = '{0}.sta{1}.{2}'.format(filebase, station, ch)
                outstream = open(fn_out, 'w')

                s1 = r"{print "
                s2 = '${0}'.format(idx_c + 1)
                s3 = r"}"
                s4 = '{0}'.format(infile)
                p = subprocess.Popen(
                    ['awk', s1 + s2 + s3, s4], stdout=outstream)
                out, dummy = p.communicate()

                outstream.close()

    metafile = open(filebase + '.timestamps', 'w')
    if samplingrate % 1 == 0:
        metafile.write('samplingrate: {0}\n'.format(int(samplingrate)))
    else:
        metafile.write('samplingrate: {0:.3f}\n'.format(samplingrate))

    if firstsample % 1 == 0:
        metafile.write('firstsample: {0}\n'.format(int(firstsample)))
    else:
        metafile.write('firstsample: {0:.3f}\n'.format(firstsample))

    if lastsample % 1 == 0:
        metafile.write('lastsample: {0}\n'.format(int(lastsample)))
    else:
        metafile.write('lastsample: {0:.3f}\n'.format(lastsample))

    metafile.close()

    if cleanup is True:
        os.remove(fileA)
        os.remove(fileB)
        os.remove(fileC)

    os.chdir(basedir)
    print('\n\t...done!\n')
    print('\t' + 72 * '-')


#========================================================

if __name__ == '__main__':
    main()
