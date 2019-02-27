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

#import pdb

lo_datafiles = ['dataA', 'dataB1', 'dataB2', 'dataC1', 'dataC2']
lo_channels = ['north', 'east', 'south', 'west', 'time']
channel_dict = {'n': 0, 'e': 1, 's': 2, 'w': 3}

scriptbase = '/wolle/elogger/test_04Jun2014/PakScriptsMay14'


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
    dataDirB = sys.argv[3]
    channelsAin = sys.argv[4]
    channelsBin = sys.argv[5]
    channelsBin = sys.argv[6]
    timestring = sys.argv[7]

    if len(sys.argv) > 8:
        dst_raw = sys.argv[8]
    try:
        dst = op.join(op.abspath(os.curdir), dst_raw)
    except:
        dst = op.join(op.abspath(os.curdir), 'birrpdata_6h')

    print()
    if not op.isdir(dst):
        os.makedirs(dst)
    for datafile in lo_datafiles:
        try:
            os.remove(op.join(dst), datafile)
        except:
            continue

    if (not op.isdir(dataDirA)):
        sys.exit("\nERROR - input data A directory not found \n")
    if (not op.isdir(dataDirB)):
        sys.exit("\nERROR - input data B directory not found \n")

    if (not op.isdir(dataDirB)):
        sys.exit("\nERROR - input data C directory not found \n")

    allchannels = 'nesw'
    try:
        channelsA = []
        tmp_chs = channelsAin.split(',')
        if len(tmp_chs) < 2:
            raise
        for i in tmp_chs:
            ch = i.lower()
            if not ch in allchannels:
                raise
            else:
                channelsA.append(channel_dict[ch])
    except:
        sys.exit("\n ERROR - unknown channels for station A\n")

    try:
        channelsB = []
        tmp_chs = channelsBin.split(',')
        if len(tmp_chs) < 2:
            raise
        for i in tmp_chs:
            ch = i.lower()
            if not ch in allchannels:
                raise
            else:
                channelsB.append(channel_dict[ch])
    except:
        sys.exit("\n ERROR - unknown channels for station B\n")
    try:
        channelsC = []
        tmp_chs = channelsCin.split(',')
        if len(tmp_chs) < 2:
            raise
        for i in tmp_chs:
            ch = i.lower()
            if not ch in allchannels:
                raise
            else:
                channelsC.append(channel_dict[ch])
    except:
        sys.exit("\n ERROR - unknown channels for station C\n")

    try:
        timestamp = time.strptime(timestring, "%y%m%d-%Hh%M")
    except:
        sys.exit("\n ERROR - unknown time - check format\n")

    # print dataDirA, dataDirB
    # print channelsA,channelsB
    # print timestamp

    # year = timestamp[0]
    # month = timestamp[1]
    # day = timestamp[2]
    # hour = timestamp[3]
    # minute=[4]

    # call the actuall processing:
    run(dataDirA, dataDirB, dataDirC, channelsA,
        channelsB, channelsC, timestamp, dst)


def run(dataDirA, dataDirB, dataDirC, channelsA,
        channelsB, channelsC, timestamp, dst):

    staAsrc, starttimeA = findDataSubDirA(timestamp, dataDirA)
    staB1src, staB2src = findDataSubDirsB(starttimeA, dataDirB)
    staC1src, staC2src = findDataSubDirsB(starttimeA, dataDirC)

    if staB1src is None and staB2src is None:
        print("\n WARNING - could not find overlapping data in given directories\n")
        pass

    if staC1src is None and staC2src is None:
        print("\n WARNING - could not find overlapping data in given directories\n")
        pass

    print(staAsrc, starttimeA, staB1src, staB2src)

    base = op.abspath(os.curdir)
    print(base)
    unpack_and_copy(staAsrc, op.join(dst, "dataA"))
    os.chdir(base)

    if staB1src is not None:
        unpack_and_copy(staB1src, op.join(dst, "dataB1"))
        os.chdir(base)
    if staB2src is not None:
        unpack_and_copy(staB2src, op.join(dst, "dataB2"))
        os.chdir(base)

    if staC1src is not None:
        unpack_and_copy(staC1src, op.join(dst, "dataC1"))
        os.chdir(base)

    if staC2src is not None:
        unpack_and_copy(staC2src, op.join(dst, "dataC2"))
        os.chdir(base)

    prepare_data(dst, channelsA, channelsB, channelsC)


def findDataSubDirsB(starttimeA, dataDirB):
    lo_subdirs = sorted(os.listdir(dataDirB))
    lo_timestamps = []
    for subdir in lo_subdirs:
        try:
            current_timestamp = calendar.timegm(
                time.strptime(subdir, "%Y-%m-%d-at-%H-%M-%S.bz"))
            lo_timestamps.append(current_timestamp)
        except:
            lo_timestamps.append(0)
            continue

    best_subdir = np.abs(starttimeA - np.array(lo_timestamps)).argmin()
    starttimeBopt = lo_timestamps[best_subdir]
    # print best_subdir, starttimeA,starttimeBopt,starttimeBopt-starttimeA

    if np.abs(starttimeBopt - starttimeA) > 3 * 3600.:
        return None, None

    if starttimeBopt < starttimeA:
        dirB1 = lo_subdirs[best_subdir]
        try:
            dirB2 = lo_subdirs[best_subdir + 1]
            starttimeB1 = starttimeBopt
            starttimeB2 = lo_timestamps[best_subdir + 1]
        except:
            dirB2 = None
        try:
            if starttimeB2 - starttimeB1 > 7 * 3600:
                raise
        except:
            dirB2 = None

    else:
        dirB2 = lo_subdirs[best_subdir]
        try:
            dirB1 = lo_subdirs[best_subdir - 1]
            starttimeB2 = starttimeBopt
            starttimeB1 = lo_timestamps[best_subdir - 1]
        except:
            dirB1 = None

        try:
            if starttimeB2 - starttimeB1 > 7 * 3600:
                raise
        except:
            dirB1 = None

    if dirB1 is not None:
        dirB1 = op.join(dataDirB, dirB1)
    if dirB2 is not None:
        dirB2 = op.join(dataDirB, dirB2)

    return dirB1, dirB2


def findDataSubDirA(timestampA, dataDirA):

    lo_subdirs = sorted(os.listdir(dataDirA))

    starttimeA = calendar.timegm(timestampA)

    lo_timestamps = []

    for subdir in lo_subdirs:
        try:
            current_timestamp = calendar.timegm(
                time.strptime(subdir, "%Y-%m-%d-at-%H-%M-%S.bz"))
            lo_timestamps.append(current_timestamp)
        except:
            lo_timestamps.append(0)
            continue

    best_subdir = np.abs(starttimeA - np.array(lo_timestamps)).argmin()

    DirA = op.join(dataDirA, lo_subdirs[best_subdir])
    starttime = lo_timestamps[best_subdir]

    return DirA, starttime


def prepare_data(dst, channelsA, channelsB, channelsC):

    base = op.abspath(os.curdir)
    os.chdir(dst)
    print('working in directory {0}'.format(os.path.abspath(os.curdir)))

    fileA = 'dataA'
    fileB1 = 'dataB1'
    fileB2 = 'dataB2'
    fileC1 = 'dataC1'
    fileC2 = 'dataC2'

    # find first and last time stamp of A
    startA = np.float64(GetFirstLine(fileA).strip().split()[-1])
    endA = np.float64(tail(fileA)[0].strip().split()[-1])

    A0_int = int(startA)
    A0_frac = int(np.round((startA - A0_int) * 1e6 / 2000)) * 2000
    A1_int = int(endA)
    A1_frac = int(np.round((endA - A1_int) * 1e6 / 2000)) * 2000

    timespanA = round(A1_int - A0_int + (A1_frac - A0_frac) * 1e-6, 5)
    no_samplesA = int(timespanA * 1e6) / 2000 + 1

    print('Info file A (from, to, samples): ', "{0:.6f}".format(startA), "{0:.6f}".format(endA), no_samplesA)
    print()

    print('writing data to single column files...')

    tA = time.gmtime(startA)
    filebase = '{0:04d}{1:02d}{2:02d}-{3:02d}h{4:02d}'.format(
        tA[0], tA[1], tA[2], tA[3], tA[4])

    print('extract time axis from file {0}'.format(fileA))
    os.system(
        r"awk '{print $5}'" +
        " {0} > {1}.time &".format(
            fileA,
            filebase))

    print('forking child 1 ...')

    pid = os.fork()
    if pid == 0:
        print('child 1 ({0}) takes over station A '.format(os.getpid()))

        for idx_c, ch in enumerate(lo_channels[:-1]):
            if idx_c in channelsA:
                print('pid ({0})'.format(os.getpid()), idx_c, ch)

                fnA = '{0}.staA.{1}'.format(filebase, ch)

                s1 = r"{print "
                s2 = '${0}'.format(idx_c + 1)
                s3 = r"}"
                s4 = '{0}'.format(fileA)
                p = subprocess.Popen(
                    ['awk', s1 + s2 + s3, s4], stdout=subprocess.PIPE)
                out, dummy = p.communicate()

                data = np.array([float(i)
                                 for i in out.split('\n') if len(i) > 0])
                with open(fnA, 'w') as F:
                    for d in data:
                        if d % 1 == 0:
                            F.write('{0}\n'.format(int(d)))
                        else:
                            F.write('{0:.4f}\n'.format(d))
                # np.savetxt(fnA,data,fmt='%.2f')
                print('file {0} done'.format(fnA))
        os._exit(0)

    if pid > 0:
        print('forking child 2 ...')
        pid2 = os.fork()
        if pid2 > 0:
            child2 = pid2
            print('parent ({0}) builds B'.format(os.getpid()))
            Bdata_section = build_section(
                fileB1,
                fileB2,
                channelsB,
                A0_int,
                A0_frac,
                A1_int,
                A1_frac,
                filebase +
                '.staB')
        else:
            print('child 2 ({0}) builds C'.format(os.getpid()))
            Cdata_section = build_section(
                fileC1,
                fileC2,
                channelsC,
                A0_int,
                A0_frac,
                A1_int,
                A1_frac,
                filebase +
                '.staC')
            os._exit(0)

        os.waitpid(child2, 0)
    for f in lo_datafiles:
        os.remove(f)

    os.chdir(base)


def build_section(fileB1, fileB2, channels, A0_int,
                  A0_frac, A1_int, A1_frac, filebase):

    startB1 = None
    endB1 = None
    values_endB1 = None

    values_startB2 = None
    startB2 = None
    endB2 = None
    timespanA = round(A1_int - A0_int + (A1_frac - A0_frac) * 1e-6, 5)
    no_samplesA = int(timespanA * 1e6) / 2000 + 1

    if op.isfile(fileB1):
        startB1 = np.float64(GetFirstLine(fileB1).strip().split()[-1])
        endB1 = np.float64(tail(fileB1)[0].strip().split()[-1])
        values_endB1 = [float(i) for i in tail(fileB1)[0].strip().split()[:-1]]
        B10_int = int(startB1)
        B10_frac = int(np.round((startB1 - B10_int) * 1e6 / 2000)) * 2000
        B11_int = int(endB1)
        B11_frac = int(np.round((endB1 - B11_int) * 1e6 / 2000)) * 2000

        timespanB1 = round(B11_int - B10_int + (B11_frac - B10_frac) * 1e-6, 5)
        no_samplesB1 = int(timespanB1 * 1e6) / 2000

    else:
        fileB1 = None

    if op.isfile(fileB2):
        values_startB2 = [float(i)
                          for i in GetFirstLine(fileB2).strip().split()[:-1]]
        startB2 = np.float64(GetFirstLine(fileB2).strip().split()[-1])
        endB2 = np.float64(tail(fileB2)[0].strip().split()[-1])
        B20_int = int(startB2)
        B20_frac = int(np.round((startB2 - B20_int) * 1e6 / 2000)) * 2000
        B21_int = int(endB2)
        B21_frac = int(np.round((endB2 - B21_int) * 1e6 / 2000)) * 2000

        timespanB2 = round(B21_int - B20_int + (B21_frac - B20_frac) * 1e-6, 5)
        no_samplesB2 = int(timespanB2 * 1e6) / 2000
    else:
        fileB2 = None

    # print startB1,endB1
    # print startB2,endB2

    if fileB1 is not None and fileB2 is not None:
        # find time axis for bridging the gap between the B files

        #...assuming 500 Hz here !!!
        # working on mu s to avoid float64 rounding issues
        t0_int = int(endB1)
        t0_frac = int(np.round((endB1 - t0_int) * 1e6 / 2000)) * 2000
        t1_int = int(startB2)
        t1_frac = int(np.round((startB2 - t1_int) * 1e6 / 2000)) * 2000

        print("end {1}: {0:.6f}".format(endB1, fileB1), "start {1}: {0:.6f}".format(startB2, fileB2))
        # print t0_int,t0_frac,t1_int,t1_frac

        gap = round(t1_int - t0_int + (t1_frac - t0_frac) * 1e-6, 5)

        # in mu s
        no_samples = int(gap * 1e6) / 2000 - 1
        # print
        # no_samples,'{0:.6f}'.format(t0_frac),'{0:.6f}'.format(t0_frac+gap*1e-6)
        taxis = (np.arange(no_samples + 1) * 2000 + t0_frac)[1:]
        taxis = taxis * 1e-6 + t0_int
        print()

        print('bridging gap  ', "{0:.6f}".format(taxis[0]), ' --> ', "{0:.6f}".format(taxis[-1]), '  ', len(taxis), 'samples')

        # and the same for the actual values...:

        bridge_data = np.zeros((len(taxis), len(values_startB2)), np.float32)
        for i in range(len(values_startB2)):
            bridge_data[
                :,
                i] = np.linspace(
                values_endB1[i],
                values_startB2[i],
                len(taxis) + 1,
                endpoint=False)[
                1:]

        # print values_endB1
        # print values_startB2
        # print bridge_data[0]
        # print bridge_data[-1]

#        sys.exit()

        # Find the origin time of A in this set of B1, B2 and bridge:
        # it must be in either B1 or the bridge
        # assume it's in B1....deal with the special case later:

        gapAB1 = round(t0_int - A0_int + (t0_frac - A0_frac) * 1e-6, 5)
        no_samplesAB1 = int(gapAB1 * 1e6) / 2000 + 1

        print('\ntake the last {0} samples from {1}'.format(no_samplesAB1, fileB1))

        gapAB2 = round(A1_int - t1_int + (A1_frac - t1_frac) * 1e-6, 5)
        no_samplesAB2 = int(gapAB2 * 1e6) / 2000 + 1

        print('take  {0} samples from gap'.format(len(bridge_data)))
        print('take the first {0} samples from {1}'.format(no_samplesAB2, fileB2))
        print('(total samples: {0})'.format(no_samplesAB1 + len(bridge_data) + no_samplesAB2))
        print()
        # sys.exit()

        dataOut = np.zeros((no_samplesA, 4), np.float32)

        #pid = os.fork()
        # if pid > 0:
        #    child = pid

        #    print 'reading file B1 ... (parent)'
        samples2skip = no_samplesB1 - no_samplesAB1 + 1
        counter = 0
        idx = 0
        for line in open(fileB1):  # in enumerate(tail(fileB1,no_samplesAB1)):
            counter += 1
            if counter <= samples2skip:
                continue
            val = line.strip().split()
            try:
                dataOut[idx] = np.array(
                    [float(val[0]), float(val[1]), float(val[2]), float(val[3])])
                if idx == 0:
                    # print val
                    pass
            except:
                print(val, idx, samples2skip)
                print(np.array([float(val[0]), float(val[1]), float(val[2]), float(val[3]), np.float64(val[4])]))
                sys.exit()

            idx += 1

        print(val)

        # sys.exit()

        # assume end is in B2
        print('interpolating gap ... ')
        dataOut[no_samplesAB1:no_samplesAB1 + len(bridge_data)] = bridge_data
        # print ' waiting for child 1'
        # os.waitpid(child, 0)

        # else:
        B2data = np.zeros((no_samplesAB2, 4), np.float32)
        print('reading file data B2/C2 ... ({0})'.format(os.getpid()))
        idx = 0
        # idx,val in enumerate(GetFirstLine(fileB2,no_samplesAB2)):
        for line in open(fileB2):
            val = line.strip().split()
            if idx == 0:
                # print val
                pass

            B2data[idx] = np.array(
                [float(val[0]), float(val[1]), float(val[2]), float(val[3])])
            idx += 1

            if idx == no_samplesAB2:
                # print val
                break
        dataOut[-no_samplesAB2:] = B2data

    for idx_c, ch in enumerate(lo_channels[:-1]):
        if idx_c in channels:
            print('pid ({0})'.format(os.getpid()), idx_c, ch)
            fnB = '{0}.{1}'.format(filebase, ch)
            data = (dataOut[:, idx_c])

            #data = detrend_linear(dataOut[:,idx_c])
            with open(fnB, 'w') as F:
                for d in data:
                    if d % 1 == 0:
                        F.write('{0}\n'.format(int(d)))
                    else:

                        F.write('{0:.10f}\n'.format(d))
                        # F.write('{0:.2f}\n'.format(d))

            print('file {0} done'.format(fnB))


def unpack_and_copy(src, dst):

    base = op.abspath(os.curdir)

    os.chdir(src)
    print(op.abspath(os.curdir))

    lo_files = sorted(os.listdir('.'))
    for f in lo_files:
        if f.lower().endswith('bz2'):
            try:
                os.system("bunzip2 -k -f {0}".format(f))
                #os.system("bzcat -k -f {0} > ".format(f))
            except:
                continue
            print('unpacked {0}'.format(f))

    pak2bin_string = op.join(scriptbase, "Pak2Bin")
    os.system('"{0}" 00000000.pak'.format(pak2bin_string))
    interp_string = op.join(scriptbase, "InterpAdl")
    os.system('"{0}" Pak2Bin-300-Data.bin o8'.format(interp_string))
    shutil.move("ELPakData.interp500.adl", dst)
    print('moving file {0} to {1}'.format("ELPakData.interp500.adl", dst))

    files = sorted(os.listdir('.'))
    files2kill = [i for i in files if i.endswith('.pak')]
    files2kill2 = [i for i in files if i.startswith('Pak2')]
    files2kill3 = [i for i in files if i.startswith('EL')]
    for f in files2kill:
        os.remove(f)
    for f in files2kill2:
        os.remove(f)
    for f in files2kill3:
        os.remove(f)

    os.chdir(base)


def GetFirstLine(fileName, n=1):
    """
    Description:        Gets last n lines using Unix tail
    Output:         returns first n lines of a file
    Keyword argument:
    n -- number of last lines to return
    filename -- Name of the file you need to tail into
    """
    p = subprocess.Popen(
        ['head', '-{0:d}'.format(n), fileName], stdout=subprocess.PIPE)
    soutput, sinput = p.communicate()

    soutput = soutput.split('\n')
    soutput = [i for i in soutput if len(i) > 0]
    if len(soutput) == 1:
        soutput = soutput[0].strip()
    # print soutput
    return soutput

    # soutput will have will contain last n lines of the code. to iterate
    # through soutput line by line do:

    # for line in GetLastNLines(50,'myfile.log').split('\n'):
    # print line


def tail(infile, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found

    lines_found = []
    f = open(infile, 'r')
    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        if len(lines_found) > lines:
            break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    output = lines_found[-lines:]
    # print output
    return output

if __name__ == '__main__':
    main()
