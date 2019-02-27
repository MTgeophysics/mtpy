#!/usr/bin/env python

"""
Fast decimation for MTpy ts-data (mtd) files
(Quick and dirty)

- no filtering applied
- only integer ratios of orignal/output sampling allowed

"""


import os
import sys
import os.path as op
import numpy as np
import mtpy.utils.filehandling as MTfh
import mtpy.utils.exceptions as MTex


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

        new_data = []
        try:
            #new_data = []
            #old_data = []
            counter = 1
            tempdata = 0
            for line in open(infile):
                line = line.strip().split()
                if line[0].startswith('#'):
                    continue
                val = float(line[0])
                tempdata += val

                if counter == 1:
                    new_data.append(val)  # tempdata/decimation_factor)
                    tempdata = 0
                counter += 1
                if counter == (decimation_factor + 1):
                    counter = 1

            # if counter != 0:
            #    new_data.append(tempdata/counter)

        # except:

                # if val%1==0:
                #    val = int(val)
                # old_data.append(val)

            #old_data = np.array(old_data)

        except:
            print('\tERROR - file does not contain single column data: {0} - SKIPPED'.format(infile))
            continue

        #len_data = len(old_data)
        if header is not None:
            n_samples = header['nsamples']

            # if len_data != n_samples:
            # print '\tWARNING - header shows wrong number of samples: {0}
            # instead of {1}'.format(n_samples,len_data)

        print('Decimating file {0} by factor {1} '.format(infile, decimation_factor))

        if n_samples % decimation_factor != 0:
            print('\tWarning - decimation of file not continuous due to mismatching decimation factor')

        #new_data = old_data[::decimation_factor]

        # index = 0
        # while index < len(new_data):
        #     old_data_index = index * decimation_factor
        #     try:
        #         new_data[index] = np.mean(old_data[old_data_index:old_data_index+decimation_factor])
        #     except:
        #         new_data[index] = np.mean(old_data[old_data_index:])

        #     index += 1

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


if __name__ == '__main__':
    run()
