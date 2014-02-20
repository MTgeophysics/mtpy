#!/usr/bin/env python

"""
Fast decimation for MTpy ts-data files
(Nasty)

- no filtering applied
- only integer ratios of orignal/output sampling allowed

"""


import os,sys
import os.path as op
import numpy as np
import mtpy.utils.filehandling as MTfh
import mtpy.utils.exceptions as MTex




def run():
    print 
    if len(sys.argv) < 4:
        sys.exit('\nNeed 3 arguments: \n\n '
            '<path to files> <output directory> <integer downsampling factor>\n \n')

    inpath = sys.argv[1]
    outpath = sys.argv[2]

    inpath = op.abspath(op.realpath(inpath))

    if not op.isdir(inpath):
        sys.exit('\nData file(s) path not existing: {0}\n'.format(inpath))

    try:
        outpath = op.abspath(op.join(os.curdir,outpath))
        if inpath == outpath:
            print 'Output directory cannot be the same as the input file location'
            raise

        if not op.exists(outpath):
            try:
                os.makedirs(outpath)
            except:
                raise
        if not os.access(outpath, os.W_OK):
            raise
    except:
        print 'Cannot generate writable output directory {0} - using'\
                ' generic location "decimated" instead'.format(outpath)
        outpath = os.path.join(inpath,'decimated')
        if not op.exists(outpath):
            os.makedirs(outpath)
        
    try:        
        decimation_factor = float(sys.argv[3])
    except:
        sys.exit('\n\tERROR - 3rd argument must be an integer decimation factor\n')

    if decimation_factor < 1 or decimation_factor%1 != 0:
        sys.exit('\n\tERROR - 3rd argument must be an integer >= 1\n')

    decimation_factor = int(decimation_factor)

    lo_files=os.listdir(inpath)
    lo_files = [i for i in lo_files if op.isfile(op.join(inpath,i))]

    if len(lo_files) == 0:
    	sys.exit('\n\tERROR - no data files in directory {0} \n'.format(inpath))
  

    for fn in lo_files:
        infile = os.path.join(inpath,fn)
        outfile = os.path.join(outpath,fn)

        try:
            header = MTfh.read_ts_header(infile)
        except MTex.MTpyError_inputarguments:
            #no TS data file
            print '\tWARNING - not a valid MTpy TS data file: {0} '.format(infile)
            continue

        old_sampling = header['samplingrate']        
        new_sampling = old_sampling/decimation_factor

        old_data = np.loadtxt(infile)
        n_samples = header['nsamples'] 

        print 'Decimating file {0} by factor {1} '.format(infile, decimation_factor)

        if n_samples%decimation_factor != 0 :
            print 'Warning - decimation of file not continuous due to mismatching decimation factor'

        new_data = old_data[::decimation_factor]

        header['nsamples'] = len(new_data)
        header['samplingrate'] = new_sampling

        new_header_line = MTfh.get_ts_header_string(header)

        Fout = open(outfile,'w')
        Fout.write(new_header_line)
        np.savetxt(Fout,new_data)
        Fout.close()

    print '\nOutput files written to {0}\n'.format(outpath)
    print '\n...Done\n'


if __name__=='__main__':
    run()

