#!/usr/bin/env python

"""
Decimation for MTpy ts-data (mtd) files


- no filtering applied
- only integer ratios of orignal/output sampling allowed

"""


import os
import sys
import os.path as op
import numpy as np
import mtpy.utils.filehandling as MTfh
import mtpy.utils.exceptions as MTex
import scipy.signal as ss


def run():
    print()
    if len(sys.argv) < 4:
        sys.exit('\nNeed 3 arguments: \n\n '
                 '<path to files> <output directory> <integer downsampling factor>\n \n')

    inpath = sys.argv[1]
    outpath = sys.argv[2]

    inpath = op.abspath(op.realpath(inpath))

    if not op.isdir(inpath):
        sys.exit('\nData file(s) path not existing: {0}\n'.format(inpath))

    try:
        outpath = op.abspath(op.join(os.curdir, outpath))
        if inpath == outpath:
            print('Output directory cannot be the same as the input file location')
            raise

        if not op.exists(outpath):
            try:
                os.makedirs(outpath)
            except:
                raise
        if not os.access(outpath, os.W_OK):
            raise
    except:
        print('Cannot generate writable output directory {0} - using'\
            ' generic location "decimated" instead'.format(outpath))
        outpath = os.path.join(inpath, 'decimated')
        if not op.exists(outpath):
            os.makedirs(outpath)

    try:
        decimation_factor = float(sys.argv[3])
    except:
        sys.exit('\n\tERROR - 3rd argument must be an integer decimation factor\n')

    if decimation_factor < 1 or decimation_factor % 1 != 0:
        sys.exit('\n\tERROR - 3rd argument must be an integer >= 1\n')

    decimation_factor = int(decimation_factor)

    lo_files = os.listdir(inpath)

    lo_files = [i for i in lo_files if op.isfile(op.join(inpath, i))]

    if len(lo_files) == 0:
        sys.exit(
            '\n\tERROR - no data files in directory {0} \n'.format(inpath))

    for fn in lo_files:
        infile = os.path.join(inpath, fn)
        outfile = os.path.join(outpath, fn)

        try:
            header = MTfh.read_ts_header(infile)
        except MTex.MTpyError_ts_data:
            # no TS data file
            print('\n\tWARNING - not a valid MTpy TS data file: {0} '.format(infile))
            header = None
            # continue

        if header is not None:
            old_sampling = header['samplingrate']
            new_sampling = old_sampling / decimation_factor

        try:
            data = []
            F = open(infile)
            for i in F:
                try:
                    if i.startswith('#'):
                        continue

                    data.append(float(i.strip()))
                except:
                    continue

            data = np.asarray(data)
            N = len(data)

        except:
            print('\tERROR - file does not contain single column data: {0} - SKIPPED'.format(infile))
            continue

        if header is not None:
            n_samples = header['nsamples']

        print('Decimating file {0} by factor {1} '.format(infile, decimation_factor))

        if n_samples % decimation_factor != 0:
            print('\tWarning - decimation of file not continuous due to mismatching decimation factor')

        # to avoid ringing in the downsampled data: use de-meaning, padding,
        # tapering:
        padlength = int(0.1 * N)
        if padlength % decimation_factor != 0:
            padlength = int(padlength / decimation_factor) * decimation_factor

        meanvalue = np.mean(data)
        data -= meanvalue
        prepad = np.ones((padlength)) * data[0]
        postpad = np.ones((padlength)) * data[-1]
        padded_data = np.append(prepad, data)
        padded_data = np.append(padded_data, postpad)

        tapered = taper_data(padded_data)

        filtered = ss.resample(tapered,
                               int(len(tapered) / decimation_factor),
                               window='blackman')

        new_padding = padlength / decimation_factor
        padding_cut = filtered[new_padding:-new_padding]

        new_data = padding_cut - np.mean(padding_cut) + meanvalue

        Fout = open(outfile, 'w')

        if header is not None:
            header['nsamples'] = len(new_data)
            header['samplingrate'] = new_sampling

            new_header_line = MTfh.get_ts_header_string(header)

            Fout.write(new_header_line)

        for i in new_data:
            if i % 1 == 0:
                Fout.write('{0:d}\n'.format(int(i)))
            else:
                Fout.write('{0:.8}\n'.format(i))
        # np.savetxt(Fout,new_data)
        Fout.close()

    print('\nOutput files written to {0}'.format(outpath))
    print('\n...Done\n')


def taper_data(data):
    """ Taper data with Tukey window function.


    input:
    -- data trace of length N

    output
    -- tapered data trace
    """

    N = len(data)

    # restriction to data with more than 50 samples
    if 0:  # (N <= 50 ):
        print(' useful tapering impossible !\n Returned original data')
        return data

    else:

        steepness = 0.85

        taper = np.ones((N), float)
        x_axis = np.arange(N) - int(N / 2)

        slopelength = int(N / 2. * (1 - 0.85) + 1)

        taper[-slopelength:] = 1. / 2. * (1 + np.cos(np.pi * (np.abs(
            x_axis[-slopelength:]) - N / 2. * steepness) / ((1 - steepness) / 2. * N)))
        taper[:slopelength] = 1. / 2. * (1 + np.cos(np.pi * (np.abs(
            x_axis[:slopelength]) - N / 2. * steepness) / ((1 - steepness) / 2. * N)))

        # for i,val in enumerate(x_axis):
        #     if N/2.*steepness <= abs(val) <= N/2.:
        #        taper[i] = 1./2.*(1+np.cos(np.pi*(np.abs(val) - N/2.*steepness)/((1-steepness)/2.*N)))

        tapered_data = data * taper

        return tapered_data  # +datamean

if __name__ == '__main__':
    run()
