#!/usr/bin/env python
"""
    MTpy Script module 

    edi2wiperdata.py


    Read EDI file(s) and extract information needed to generate a PhaseTensorWiper plot 
    representation of data

    Actual PT Wiper Plot is then generated using GMT
    (to be included later)


    @MTpy 2014, UofA (LK)

"""


import os,sys
import os.path as op

import mtpy.core.edi as MTedi
import mtpy.analysis.geometry as MTgy
#import mtpy.analysis.pt as MTpt
import mtpy
#for debugging:

import pdb
#reload(MTedi)
#reload(MTpt)
reload(MTgy)


def main():


    if len(sys.argv) < 3:
        print '\nNeed at least 2 arguments: <EDI file> '\
                        '<output directory> \n\n'\
                        'Optional arguments: \n [output filename]\n'\
                        ' [batch process flag "-b"] \n\n'
        return

    try:
        fn_in = sys.argv[1]
        fn_in = op.join(op.abspath(os.curdir),fn_in)
        edi_object = MTedi.Edi(filename=fn_in)
    except:
        print '\n\tERROR - File is not a valid EDI file: {0}\n'.format(fn_in)
        sys.exit()

    try:
        outdir = sys.argv[2]
        outdir = op.join(op.abspath(os.curdir),outdir)
        if not op.isdir(outdir):
            os.makedirs(outdir)
    except:
      print '\n\tERROR - Output directory does not exist and cannot be'\
                                          ' generated: {0}\n'.format(outdir)
      sys.exit()

    fn_out = None
    if len(sys.argv)>3:
        try:
            fn_out = sys.argv[3]
            fn_out = op.join(outdir,fn_out)
        except:
            fn_out = None
            print '\nWARNING - Invalid output filename...using default instead\n'

    try:
        generate_ptwiperdata_file(edi_object, outdir, fn_out)
    except:
        print '\n\tERROR - could not generate PT Wiper Data file - check EDI file!\n'


def generate_ptwiperdata_file(edi_object, outdir,outfn=None):

    pt = mtpy.analysis.pt.PhaseTensor(z_object=edi_object.Z,freq=edi_object.freq)

    #debugging stop:
    #pdb.set_trace()

    station = edi_object.station
    #no spaces in file names:
    if len(station.split())>1:
        station = '_'.join(station.split())

    a = pt.alpha
    b = pt.beta
    pmin = pt.phimin
    pmax = pt.phimax
    e = pt.ellipticity

    if outfn is None:
        fn = '{0}_PTwiperdata'.format(station)
        outfn = op.join(outdir,fn)
    
    outfn = op.realpath(outfn)

    try:
        Fout = open(outfn,'w')
    except:
        print '\n\tERROR - Cannot generate output file!\n'
        raise
    Fout.write('# {0}   {1:+010.6f}   {2:+011.6f}\n'.format(station,edi_object.lat,edi_object.lon))
    headerstring = '# freq \t\t Phimin  sigma \t Phimax  sigma \t alpha  '\
                    'sigma \t beta  sigma \t ellipticity  sigma \n'
    Fout.write(headerstring)
    for i,freq in enumerate(edi_object.freq):
        try:
            vals = '{0:.4e}\t{1: 3.2f}\t{2:3.2f}\t{3: 3.2f}\t{4:3.2f}\t{5: 3.2f}\t{6:3.2f}'\
            '\t{7: 3.2f}\t{8:3.2f}\t{9:.3f}\t{10:.3f}\n'.format(
                freq,pmin[0][i],pmin[1][i],pmax[0][i],pmax[1][i],a[0][i],a[1][i],
                b[0][i],b[1][i],e[0][i],e[1][i])            
            Fout.write(vals)
        except:
            continue
    
    Fout.close()
    print '\n\t Done - Written data to file: {0}\n'.format(outfn)

if __name__=='__main__':
    main()


