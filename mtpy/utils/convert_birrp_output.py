#!/usr/bin/env python
"""
mtpy/mtpy/uofa/convert_birrp_output.py

This is a convenience script for converting BIRRP output data into EDI- and 
coherence-files. 

arguments:

stationname, input directory, survey_configfile, birrp_configfile

optional:
instrument response file  


"""



import sys
import os
import os.path as op


import mtpy.processing.birrp as MTbp
reload(MTbp)


def main():

    if len(sys.argv) < 5:
        sys.exit('\nNeed at least 4 arguments:\n '
                '<stationname> \n <path to Birrp output files> \n '
                '<survey config file>\n <Birrp config file>\n\n'
                '[optional: instrument response file in (freq,real,imag)-format]\n')
        
    stationname = sys.argv[1]
    datadir = sys.argv[2]
    survey_cfg_fn = sys.argv[3]
    birrp_cfg_fn = sys.argv[4]
    
    instr_resp_fn = None
    if len(sys.argv) > 5:
        instr_resp_fn = sys.argv[5]


    convertbirrp(stationname,datadir,survey_cfg_fn,birrp_cfg_fn,instr_resp_fn)

    

def convertbirrp(stationname, datadir, survey_configfile, 
                birrp_configfile, instr_response_file):


    try:
        datadir = op.realpath(op.abspath(datadir))
        if not op.isdir(datadir):
            raise
    except:
        sys.exit('Directory not existing: {0}'.format(datadir))

    try:
        survey_configfile = op.abspath(survey_configfile)
        if not op.isfile(survey_configfile):
            raise
    except:
        sys.exit('Survey config file not existing: {0}'.format(survey_configfile))

    try:
        birrp_configfile = op.abspath(birrp_configfile)
        if not op.isfile(birrp_configfile):
            raise
    except:
        sys.exit('Birrp config file not existing: {0}'.format(birrp_configfile))

    if instr_response_file is not None:
        try:
            ir_fn = op.abspath(instr_response_file)
            if not op.isfile(ir_fn):
                raise
        except:
            sys.exit('Instrument response file not existing: {0}'.format(ir_fn))

        try:
            MTbp.convert2coh(stationname, datadir)
        except:
            try:
                MTbp.convert2coh(stationname.upper(), datadir) 
            except:
                print 'Could not generate coherence file'

        try:
            MTbp.convert2edi_incl_instrument_correction(stationname,\
                                                        datadir,\
                                                        survey_configfile,\
                                                        birrp_configfile,\
                                                        ir_fn)
        except:
            try:
                MTbp.convert2edi_incl_instrument_correction(stationname.upper(),\
                                                        datadir,\
                                                        survey_configfile,\
                                                        birrp_configfile,\
                                                        ir_fn)
            except:
                print 'Could not generate EDI file'
        
        return
 

    try:
        MTbp.convert2coh(stationname, datadir)
    except:
        try:
            MTbp.convert2coh(stationname.upper(), datadir) 
        except:
            print 'Could not generate coherence file'

    try:
        MTbp.convert2edi(stationname,\
                        datadir,\
                        survey_configfile,\
                        birrp_configfile)
    except:
        try:
            MTbp.convert2edi(stationname.upper(),\
                            datadir,\
                            survey_configfile,\
                            birrp_configfile)
        
        except:
            print 'Could not generate EDI file'
    
    

if __name__=='__main__':
    main()
