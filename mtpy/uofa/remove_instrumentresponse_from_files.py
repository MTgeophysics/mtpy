#!/usr/bin/env python
"""
This is a convenience script for the removal of instrument response from a set of MTpy time series data files within a directory (non-recursive). The data files have to contain a MTpy style header line, which specifies station, channel, timestamps.

It needs the location of the directory and the location of the instrument response file.
The latter has to consist of an array with three columns: frequencies, real, imaginary

If no output folder is specified, a subfolder 'instr_resp_corrected' is set up within the input directory

"""

import numpy as np
import re
import sys
import os
import glob
import os.path as op
import glob
import calendar
import time
import scipy.signal


import mtpy.utils.exceptions as MTex
import mtpy.processing.calibration as MTcb
import mtpy.utils.filehandling as MTfh
import mtpy.processing.instrument as MTin

#reload(MTfh)
#reload(MTex)
#reload(MTcb)
#reload(MTin)


def main():

    if len(sys.argv) < 3:
        raise MTex.MTpyError_inputarguments(
            'Need at least 2 arguments: <path to files> <response file> [<output dir>] [<channel(s)>] ')

    pathname_raw = sys.argv[1]
    directory = op.abspath(op.realpath(pathname_raw))

    responsefilename_raw = sys.argv[2]
    responsefile = op.abspath(op.realpath(responsefilename_raw))

    if not op.isdir(directory):
        raise MTex.MTpyError_inputarguments(
            'Directory not existing: %s' % (directory))

    if not op.isfile(responsefile):
        raise MTex.MTpyError_inputarguments(
            'Response file not existing: %s' %
            (responsefile))

    # check, if response file is in proper shape (3 columns freq,re,im of real
    # values):
    try:
        responsedata = np.loadtxt(responsefile)
        s = responsedata.shape
        if s[1] != 3:
            raise
        freq_min = responsedata[0, 0]
        freq_max = responsedata[-1, 0]

    except:
        raise MTex.MTpyError_inputarguments(
            'Response file (%s) in wrong format - must be 3 columns: freq,real,imag' %
            (responsefile))

    # set up output directory:
    try:
        outdir_raw = sys.argv[3]
        outdir = op.abspath(outdir_raw)
    except:
        outdir = op.join(directory, 'instr_resp_corrected')

    try:
        if not op.isdir(outdir):
            os.makedirs(outdir)
    except:
        raise MTex.MTpyError_inputarguments(
            'Output directory cannot be generated: %s' %
            (outdir))

    # define channels to be considered for correction:
    try:
        lo_channels = list(set([i.upper() if len(
            i) == 2 else 'B' + i.upper() for i in sys.argv[4].split(',')]))
    except:
        print('No channel list found - using BX, BY, HX, HY')
        lo_channels = ['BX', 'BY', 'HX', 'HY', 'BZ', 'HZ']

    # collect file names  within the folder
    lo_allfiles = os.listdir(directory)

    lo_allfiles = [
        op.abspath(
            op.join(
                directory,
                i)) for i in lo_allfiles if op.isfile(
            op.abspath(
                op.join(
                    directory,
                    i)))]  # if op.isfile(i)==True]

    # generate list of list-of-files-for-each-channel:
    lo_lo_files_for_channels = [[] for i in lo_channels]

    # check the files for information about the determined channels:
    for fn in lo_allfiles:
        header_dict = MTfh.read_ts_header(fn)
        if len(list(header_dict.keys())) == 0:
            continue
        ch = header_dict['channel'].upper()
        if ch not in lo_channels:
            continue

        try:
            ch_idx = lo_channels.index(ch)
        except ValueError:
            continue

        # use the current file, if it contains a header line and contains
        # signal from the requested channel:
        lo_lo_files_for_channels[ch_idx].append(fn)

    # if no files had header lines or did not contain data from the
    # appropriate channel(s):
    if np.sum([len(i) for i in lo_lo_files_for_channels]) == 0:
        print('channels: ', lo_channels, ' - directory: ', directory)
        raise MTex.MTpyError_inputarguments(
            'No information for channels found in directory {0} - Check header lines!'.format(directory))

    #=============================================
    # start the instrument correction

    # looping over all requested channels:
    for ch in lo_channels:
        # skip, if no data are available for the current channel:
        if [len(i) for i in lo_lo_files_for_channels][
                lo_channels.index(ch)] == 0:
            continue

        # set up lists for the infos needed later, esp. for the file handling
        lo_files = lo_lo_files_for_channels[lo_channels.index(ch)]
        lo_t_mins = []
        lo_headers = []

        # read in header lines and sort files by increasing starttimes t_min
        for fn in lo_files:
            header_dict = MTfh.read_ts_header(fn)
            lo_t_mins.append(header_dict['t_min'])
            lo_headers.append(header_dict)

        # sort all the collected lists by t_min
        idxs = np.array(lo_t_mins).argsort()
        lo_t_mins = [lo_t_mins[i] for i in idxs]
        lo_files = [lo_files[i] for i in idxs]
        lo_headers = [lo_headers[i] for i in idxs]

        # finding consecutive, continuous time axes:
        lo_timeaxes = []
        ta_old = None

        for idx, header in enumerate(lo_headers):
            ta_cur = np.arange(
                int(header['nsamples'])) / float(header['samplingrate']) + float(header['t_min'])

            # if there is no old ta:
            if ta_old is None:
                ta_old = ta_cur
                continue

            # if gap between old and new ta is too big:
            if (ta_cur[0] - ta_old[-1]) > (2 *
                                           1. / float(header['samplingrate'])):
                lo_timeaxes.append(np.array(ta_old))
                ta_old = ta_cur
                continue

            # find index of new ta which is closest to the end of old_ta - most
            # commonly it's '0' !
            overlap = np.abs(ta_cur - ta_old[-1]).argmin()
            ta_cur = ta_cur[overlap:]
            ta_old = np.concatenate([ta_old, ta_cur])

        # append last active time axis ta:
        lo_timeaxes.append(np.array(ta_old))

        # determine maximal period from response file and existinng time axes.
        # win = get_windowlength() = max([ i[-1]-i[0] for i in lo_timeaxes] )
        # the minimum of the maximal resolvable signal period and the longest
        # continuous time axis:
        winmax = 1. / freq_min
        # for debugging set large window size:
        #winmax = 5e5
        # later on, if the TS is longer than 3 times this time window, we want
        # to cut out subsections of the time series. These cuts shall consist
        # of triplets of subwindows, each of which shall not be longer than
        # this maximum period.

        # Now the data set has to be corrected/deconvolved by looping over the
        # collected time axes:
        for ta in lo_timeaxes:
            print('\nhandling time axis: {0} - {1} ({2} samples) '.format(ta[0], ta[-1], len(ta)))

            # if the time section is shorter than 3 times the maximum defined
            # by the response function, read in the whole data set at once for
            # this interval
            if (ta[-1] - ta[0]) < (3 * winmax):
                print('time axis short enough ({0} seconds) - reading all at once'.format(ta[-1] - ta[0]))

                # collect the appropriate files in a list
                # after the MTpy preprocessing the start end end of the time
                # series coincide with files start and endings, so no "half
                # files" are involved.
                cur_time = ta[0]
                data = []
                files = []
                headers = []
                starttimes = []

                while cur_time < ta[-1]:
                    for idx, header in enumerate(lo_headers):
                        ta_cur = np.arange(
                            int(header['nsamples'])) / float(header['samplingrate']) + float(header['t_min'])
                        if cur_time in ta_cur:
                            start_idx = np.where(ta_cur == cur_time)[0][0]
                            break
                    fn = lo_files[idx]
                    files.append(fn)
                    headers.append(header)
                    starttimes.append(float(header['t_min']))
                    cur_data = np.loadtxt(fn)

                    print('current data section length: ', len(cur_data))
                    if ta_cur[-1] <= ta[-1]:
                        data.extend(cur_data[start_idx:].tolist())
                        cur_time = ta_cur[-1] + 1. / \
                            float(header['samplingrate'])
                    else:
                        end_idx = np.where(ta_cur == ta[-1])[0][0]
                        data.extend(cur_data[start_idx:end_idx + 1].tolist())
                        cur_time = ta[-1]
                    print('current data length: ', len(data))

                data = np.array(data)
                data = scipy.signal.detrend(data)
                # at this point, the data set should be set up for the given
                # time axis
                corrected_timeseries = MTin.correct_for_instrument_response(
                    data, float(header['samplingrate']), responsedata)

                print('corrected TS starting at {0}, length {1}'.format(ta[0], len(corrected_timeseries)))

                # now, save this TS back into the appropriate files, including
                # headers
                for idx, fn in enumerate(files):

                    # output file name: use input file name and append '_true'
                    inbasename = op.basename(fn)
                    outbasename = ''.join(
                        [op.splitext(inbasename)[0] + '_true', op.splitext(inbasename)[1]])
                    outfn = op.join(outdir, outbasename)

                    outF = open(outfn, 'w')
                    header = headers[idx]
                    unit = header['unit']
                    if unit[-6:].lower() != '(true)':
                        unit += '(true)'
                    header['unit'] = unit
                    headerline = MTfh.get_ts_header_string(header)
                    outF.write(headerline)
                    starttime = starttimes[idx]
                    length = int(float(header['nsamples']))

                    startidx = (np.abs(starttime - ta)).argmin()
                    print(startidx, length, len(corrected_timeseries), len(ta))
                    print('\nhandling file {0} - starttime {1}, - nsamples {2}'.format(outfn, starttime, length))
                    print(outdir, outfn)
                    data = corrected_timeseries[startidx:startidx + length]
                    np.savetxt(outF, data)
                    outF.close()

                # To do so, use the time axis and run over the input files again,determine the filenames. Use them, and put them (slightly modified, perhaps?) into the given output directory
                # return corrected_timeseries

            else:

                # find partition into pieces of length 'winmax'. the remainder
                # is equally split between start and end:

                # assume constant sampling rate, just use the last opened
                # header (see above):
                samplingrate = float(header['samplingrate'])

                # total time axis length:
                ta_length = ta[-1] - ta[0] + 1. / samplingrate

                # partition into winmax long windows
                n_windows = int(ta_length / winmax)
                remainder = ta_length % winmax
                lo_windowstarts = [ta[0]]
                for i in range(n_windows + 1):
                    t0 = ta[0] + remainder / 2. + i * winmax
                    lo_windowstarts.append(t0)
                print('time axis long ({0} seconds) - processing in {1} sections (window: {2})'.format(ta_length, n_windows, winmax))

                # lists of input file(s) containing the data - for all 3
                # sections of the moving window
                section1_lo_input_files = []
                section2_lo_input_files = []
                section3_lo_input_files = []
                section1_data = []
                section2_data = []
                section3_data = []
                section1_ta = []
                section2_ta = []
                section3_ta = []
                file_open = False

                # loop over the winmax long sections:
                for idx_t0, t0 in enumerate(lo_windowstarts):

                    # the last two window-starts do not get an own window:
                    if idx_t0 > n_windows - 1:
                        # since there are 'n_windows' full intervals for the
                        # moving window
                        break

                    print('\n----- section {0} ----\n'.format(idx_t0 + 1))
                    if idx_t0 == n_windows - 1:
                        endtime = ta[-1]
                    else:
                        endtime = lo_windowstarts[
                            idx_t0 + 3] - 1. / samplingrate
                    print('s1 = {0} - {1}, s2 = {2} - {3}, s3 = {4} - {5}\n'.format(t0, lo_windowstarts[idx_t0 + 1] - 1. / samplingrate, lo_windowstarts[idx_t0 + 1], lo_windowstarts[idx_t0 + 2] - 1. / samplingrate, lo_windowstarts[idx_t0 + 2], endtime))

                    # for each step (except for he last one obviously), 3 consecutive parts are read in, concatenated and deconvolved. Then the central part is taken as 'true' data.
                    # only for the first and the last sections (start and end
                    # pieces) are handled together with the respective
                    # following/preceding section

                    # the data currently under processing:
                    data = []
                    timeaxis = []

                    # list of current input data files and their starting times
                    # - only files, which contain parts of the time interval
                    # for the output!
                    lo_infiles = []
                    lo_t0s = []

                    # if old data are present from the step before:
                    if (len(section2_data) > 0) and (len(section3_data) > 0):
                        print('old data found....moving window')
                        section1_data = section2_data
                        section2_data = section3_data
                        section1_ta = section2_ta
                        section2_ta = section3_ta
                        section1_lo_input_files = section2_lo_input_files
                        section2_lo_input_files = section3_lo_input_files
                        section1_lo_t0s = section2_lo_t0s
                        section2_lo_t0s = section3_lo_t0s

                    # otherwise, it's the first step, so all 3 sections have to
                    # be read
                    else:
                        print('first section...initialise data collection')
                        section1_data, section1_lo_input_files, section1_lo_t0s = read_ts_data_from_files(
                            t0, lo_windowstarts[idx_t0 + 1], lo_t_mins, lo_files)
                        section2_data, section2_lo_input_files, section2_lo_t0s = read_ts_data_from_files(
                            lo_windowstarts[idx_t0 + 1], lo_windowstarts[idx_t0 + 2], lo_t_mins, lo_files)

                        section1_ta = np.arange(
                            len(section1_data)) / samplingrate + t0
                        section2_ta = np.arange(
                            len(section2_data)) / samplingrate + lo_windowstarts[idx_t0 + 1]

                    if idx_t0 < n_windows - 1:
                        print('lll')
                        section3_data, section3_lo_input_files, section3_lo_t0s = read_ts_data_from_files(
                            lo_windowstarts[idx_t0 + 2], lo_windowstarts[idx_t0 + 3], lo_t_mins, lo_files)
                    else:
                        print('jjjj')
                        # for the last section, there is no
                        # lo_windowstarts[idx_t0 +3], so it must be the end of
                        # the overall time axis
                        section3_data, section3_lo_input_files, section3_lo_t0s = read_ts_data_from_files(
                            lo_windowstarts[idx_t0 + 2], ta[-1] + 1. / samplingrate, lo_t_mins, lo_files)

                    section3_ta = np.arange(
                        len(section3_data)) / samplingrate + lo_windowstarts[idx_t0 + 2]
                    data = np.concatenate(
                        [section1_data, section2_data, section3_data])
                    timeaxis = np.concatenate(
                        [section1_ta, section2_ta, section3_ta])
                    print('sections ta: ', section1_ta[0], section1_ta[-1], len(section1_ta), section2_ta[0], section2_ta[-1], len(section2_ta), section3_ta[0], section3_ta[-1], len(section3_ta))
                    # continue
                    data = np.array(data)
                    # remove mean and linear trend (precaution, since it should be removed as low frequency
                    # content anyway)
                    data = scipy.signal.detrend(data)

                    #----------------------------
                    # the actual deconvolution:
                    corrected_data = MTin.correct_for_instrument_response(
                        data, samplingrate, responsedata)
                    #----------------------------

                    if idx_t0 == 0:
                        # for the first window, the first section is output as
                        # well
                        startidx = 0
                        lo_infiles.extend(section1_lo_input_files)
                        lo_t0s.extend(section1_lo_t0s)

                    else:
                        # otherwise, just the central section, so take the
                        # start of this middle section
                        startidx = (
                            np.abs(
                                timeaxis -
                                lo_windowstarts[
                                    idx_t0 +
                                    1])).argmin()

                    # collect the respective input filenames and their starting
                    # times
                    for fn in section2_lo_input_files:
                        if fn not in lo_infiles:
                            lo_infiles.append(fn)
                    for t0 in section2_lo_t0s:
                        if t0 not in lo_t0s:
                            lo_t0s.append(t0)

                    if idx_t0 == n_windows - 1:
                        # for the last window, the last section is added
                        for fn in section3_lo_input_files:
                            if fn not in lo_infiles:
                                lo_infiles.append(fn)
                        for t0 in section3_lo_t0s:
                            if t0 not in lo_t0s:
                                lo_t0s.append(t0)

                        data2write = corrected_data[startidx:]
                        timeaxis2write = timeaxis[startidx:]

                    else:
                        # for all windows but the last, get the middle section
                        # end
                        endidx = (
                            np.abs(
                                timeaxis -
                                lo_windowstarts[
                                    idx_t0 +
                                    2])).argmin()
                        # and cut out the section
                        data2write = corrected_data[startidx:endidx]
                        timeaxis2write = timeaxis[startidx:endidx]

                    # print
                    # 'indizes:',len(timeaxis),startidx,len(timeaxis2write),timeaxis2write[0],timeaxis2write[-1]

                    # at this point, we got
                    # - the corrected data
                    # - the sorted list of the involved input files (and their starting times)
                    #
                    # write data to the output files, whose structure is
                    # following the input file(s):

                    # maximum time stamp of the current data set (section)
                    tmax = timeaxis2write[-1]

                    # initialise the current time to the beginning of the
                    # section to be written
                    t = timeaxis2write[0]
                    print('\nfiles involved in this section', lo_infiles)

                    while t < tmax:
                        print(t)

                        if file_open == False:
                            print('no file open...preparing new one')
                            # take the first of the input files in the list:
                            print('read header of input file {0} '.format(lo_infiles[0]))
                            header = MTfh.read_ts_header(lo_infiles[0])
                            ta_tmp = np.arange(float(
                                header['nsamples'])) / float(header['samplingrate']) + float(header['t_min'])
                            unit = header['unit']
                            if unit[-6:].lower() != '(true)':
                                unit += '(true)'
                            header['unit'] = unit
                            headerline = MTfh.get_ts_header_string(header)

                            # output file name: use input file name and append
                            # '_true'
                            inbasename = op.basename(fn)
                            outbasename = ''.join(
                                [op.splitext(inbasename)[0] + '_true', op.splitext(inbasename)[1]])
                            outfn = op.join(outdir, outbasename)
                            print('write header to output file {0} '.format(outfn))
                            outF = open(outfn, 'w')
                            outF.write(headerline)

                            # if the section exceeds the time axis of the file:
                            if tmax > ta_tmp[-1]:
                                print('data longer than space left in file...storing first part ')
                                # write as many samples to the files as there
                                # belong
                                np.savetxt(
                                    outF, data2write[
                                        :int(
                                            float(
                                                header['nsamples']))])
                                print('write data to output file {0} '.format(outfn))

                                # close the file
                                outF.close()
                                file_open = False
                                print(' output file closed : {0} '.format(outfn))
                                datalength = 0
                                # drop out the first elements of the lists
                                dummy = lo_infiles.pop(0)
                                print('dropped {0}'.format(dummy))
                                dummy = lo_t0s.pop(0)
                                print('dropped {0}'.format(dummy))

                                # cut the written part of the data
                                print('data cut from {0} to {1}'.format(len(data2write), len(data2write[int(float(header['nsamples'])):])))
                                data2write = data2write[
                                    int(float(header['nsamples'])):]
                                timeaxis2write = timeaxis2write[
                                    int(float(header['nsamples'])):]
                                # define the current time as one sample after
                                # the end of the file, which was just closed
                                t = ta_tmp[-1] + 1. / \
                                    float(header['samplingrate'])
                                print('current time set to {0}'.format(t))
                                # and back to the while condition, since there
                                # are unwritten data

                            # if the section is not longer than the time axis
                            # of the newly opened file:
                            else:
                                print('data fits into open file : {0}'.format(outF.name))
                                # write everything
                                np.savetxt(outF, data2write)
                                print('wrote data to file')
                                # check, if by chance this is exactly the
                                # correct number of samples for this file:
                                if tmax == ta_tmp[-1]:
                                    # if so, close it
                                    print('file full....closing...')
                                    outF.close()
                                    file_open = False
                                    datalength = 0

                                    print('written data to file and closed it: {0}'.format(outF.name))

                                    dummy = lo_infiles.pop(0)
                                    dummy = lo_t0s.pop(0)
                                    # actually, the infile list should be empty
                                    # now, since the whole section has been
                                    # written to a file!!

                                # otherwise, the section is shorter, so the
                                # file (potentially) misses entries
                                else:
                                    datalength = len(data2write)

                                    print('file not full...{0} remains open - {1} of max. {2} samples written '.format(outF.name, datalength, int(float(header['nsamples']))))
                                    file_open = True

                                # define the current time as at the end of the
                                # time axis of the section, i.e. 'go to next
                                # section':
                                t = tmax
                                print('current time set to ', t)

                        # otherwise, a file is already open and the next
                        # section has to be appended there:
                        else:
                            print('open file {0} is waiting for data...'.format(outF.name))
                            header = MTfh.read_ts_header(lo_infiles[0])
                            ta_tmp = np.arange(float(
                                header['nsamples'])) / float(header['samplingrate']) + float(header['t_min'])
                            print(tmax, ta_tmp[-1], '{0} of {1} samples used...{2} waiting\n'.format(datalength, len(ta_tmp), len(data2write)))

                            # check, if the data exceeds the time axis of the
                            # open file:
                            if tmax > ta_tmp[-1]:
                                print('data longer than space left in file...storing first {0} here and cut of the rest ({1}) for further processing...\n'.format(len(ta_tmp) - datalength, len(data2write) - len(ta_tmp) + datalength))

                                # determine the index of the section time axis,
                                # which belongs to the last entry of the
                                # currently open file - including the last
                                # value!
                                endidx = (np.abs(timeaxis2write -
                                                 ta_tmp[-1])).argmin() + 1
                                # write the respective part of the data to the
                                # file and close it then
                                np.savetxt(outF, data2write[:endidx])
                                outF.close()
                                file_open = False
                                print('closed full file {0}'.format(outF.name))
                                # cut out the first bit, which is already
                                # written
                                print('cutting data to remaining bit ({0}->{1})\n'.format(len(data2write), len(data2write) - len(data2write[:endidx])))
                                data2write = data2write[endidx:]
                                timeaxis2write = timeaxis2write[endidx:]
                                # drop the file, which is used and done
                                print('dropping name of closed file', lo_infiles[0])
                                dummy = lo_infiles.pop(0)
                                dummy = lo_t0s.pop(0)
                                # set the current time to the start of the next
                                # file
                                t = ta_tmp[-1] + 1. / \
                                    float(header['samplingrate'])
                                print('current time set to ', t)

                            # if the section is not longer than the time axis
                            # of the open file:
                            else:
                                print('data fits into the open file :', outF.name)
                                # write everything
                                np.savetxt(outF, data2write)
                                datalength += len(data2write)
                                print('file contains {0} samples - {1} open '.format(datalength, int(float(header['nsamples'])) - datalength))
                                # check, if by chance this is exactly the
                                # correct number of samples for this file:
                                if tmax == ta_tmp[-1]:
                                    # if so, close it
                                    outF.close()
                                    print('written data to file {0}'.format(outfn))
                                    file_open = False
                                    dummy = lo_infiles.pop(0)
                                    dummy = lo_t0s.pop(0)
                                    # actually, the infile list should be empty
                                    # now, since the whole section has been
                                    # written to a file!!

                                # otherwise, the section is shorter, so the
                                # file (potentially) misses entries
                                else:
                                    file_open = True

                                # define the current time as at the end of the
                                # time axis of the section => go to next
                                # section:
                                t = tmax

                # after looping over all sections, check, if the last file has
                # been closed:
                if file_open == True:
                    outF.close()
                    print('written data to file {0}'.format(outF.name))


def read_ts_data_from_files(
        starttime, endtime, list_of_starttimes, list_of_files):
    """
        Use the tmin_list to determine, which files have to be opened to obtain data from the time t0 (inclusive) to 'endtime(exclusive!!)'

        Return the data as well as a list of files
    """

    # sorting lists for starttimes:
    sorted_idx = np.array(list_of_starttimes).argsort()
    list_of_starttimes = [float(list_of_starttimes[i]) for i in sorted_idx]
    list_of_files = [list_of_files[i] for i in sorted_idx]
    filelist = []
    t0_list = []
    data = []

    if starttime < list_of_starttimes[0]:
        sys.exit(
            'Cannot read any data - requested interval is not covered by the data files time axis -starttime {0} - {1}'.format(
                starttime,
                list_of_starttimes))

    file_idx = 0
    t = starttime
    while t < endtime:
        fn = list_of_files[file_idx]
        header = MTfh.read_ts_header(fn)
        ta_tmp = np.arange(float(
            header['nsamples'])) / float(header['samplingrate']) + float(header['t_min'])
        if ta_tmp[0] <= t <= ta_tmp[-1]:
            startidx = (np.abs(t - ta_tmp)).argmin()

            if ta_tmp[0] <= endtime <= ta_tmp[-1]:
                endidx = (np.abs(endtime - ta_tmp)).argmin()
                data.extend(list(np.loadtxt(fn)[startidx:endidx]))
                t = endtime
            else:
                data.extend(list(np.loadtxt(fn)[startidx:]))
                t = ta_tmp[-1] + 1. / float(header['samplingrate'])
            if fn not in filelist:
                filelist.append(fn)
            t0_list.append(float(header['t_min']))
        file_idx += 1

    return np.array(data), filelist, t0_list


if __name__ == '__main__':
    main()
