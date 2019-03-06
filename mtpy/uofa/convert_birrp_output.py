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
import numpy as np

import mtpy.processing.birrp as MTbp
#reload(MTbp)


def main():

    if len(sys.argv) < 4:
        sys.exit('\nNeed at least 3 arguments:\n '
                 '<stationname> \n <path to Birrp output files> \n '
                 '<survey config file>\n\n[optinal: -b <Birrp config file>\n'
                 '[optional: -i <instrument response file> in (freq,real,imag)-format]\n')

    stationname = sys.argv[1]
    datadir = sys.argv[2]
    survey_cfg_fn = sys.argv[3]

    birrp_cfg_fn = None
    instr_resp_fn = None

    if len(sys.argv) > 4:
        optionals = sys.argv[4:]

        for idx_o, o in enumerate(optionals):
            if o[0] == '-':
                option = o[1].lower()
                if option not in ['b', 'i']:
                    print('unknown option: {0}'.format(option))
                    continue
                else:
                    try:
                        argument = ''
                        argument = optionals[idx_o + 1]
                        if argument[0] == '-':
                            argument = ''
                            raise
                        if option == 'b':
                            birrp_cfg_fn = argument
                        if option == 'i':
                            instr_resp_fn = argument
                    except:
                        print('option "{0}" not followed by valid argument: "{1}"'\
                            ''.format(option, argument))

    edifn, cohfn = convertbirrpoutput(
        stationname, datadir, survey_cfg_fn, birrp_cfg_fn, instr_resp_fn)

    print('EDI/coh - files generated for station {0}:\n{1}\n{2}'\
        ''.format(stationname, edifn, cohfn))


def convertbirrpoutput(stationname, datadir, survey_configfile, birrp_configfile=None,
                       instr_response_file=None):

    edifn = None
    cohfn = None

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
        sys.exit(
            'Survey config file does not exist: {0}'.format(survey_configfile))

    if birrp_configfile is not None:
        try:
            birrp_configfile = op.abspath(birrp_configfile)
            if not op.isfile(birrp_configfile):
                raise
        except:
            print('Birrp config file not existing: {0} - using generic values'\
                ''.format(birrp_configfile))
            birrp_configfile = None

    if instr_response_file is not None:
        try:

            ir_fn = op.abspath(instr_response_file)
            # print irfn,instr_response_file
            if not op.isfile(ir_fn):
                raise
        except:
            sys.exit(
                'Instrument response file does not exist: {0}'.format(ir_fn))
        try:
            instr_resp = np.loadtxt(ir_fn)
            if np.shape(instr_resp)[1] != 3:
                raise
            if len(np.shape(instr_resp)) != 2:
                raise
            if np.shape(instr_resp)[0] < 2:
                raise
        except:
            sys.exit('\n\t!!! Instrument response file has wrong format !!!\n'
                     '\nNeeds 3 columns and at least 2 rows containing complex'
                     ' valued transfer function:\n\n\tfrequency, real,'
                     ' imaginary\n\n')

        try:
            cohfn = MTbp.convert2coh(stationname, datadir)
        except:
            try:
                print('trying to find files for uppercase stationname')
                cohfn = MTbp.convert2coh(stationname.upper(), datadir)
            except:
                cohfn = None
                print('Could not generate coherence file')
        try:
            edifn = MTbp.convert2edi_incl_instrument_correction(stationname,
                                                                datadir,
                                                                survey_configfile,
                                                                birrp_configfile,
                                                                ir_fn)

        except:
            raise
            try:
                edifn = MTbp.convert2edi_incl_instrument_correction(stationname.upper(),
                                                                    datadir,
                                                                    survey_configfile,
                                                                    birrp_configfile,
                                                                    ir_fn)
            except:
                edifn = None
                print('Could not generate EDI file')

        return edifn, cohfn

    try:
        cohfn = MTbp.convert2coh(stationname, datadir)
    except:
        try:
            print('trying to find files for uppercase stationname')
            cohfn = MTbp.convert2coh(stationname.upper(), datadir)
        except:
            cohfn = None
            print('Could not generate coherence file')

    try:
        edifn = MTbp.convert2edi(stationname,
                                 datadir,
                                 survey_configfile,
                                 birrp_configfile)
    except:
        try:
            edifn = MTbp.convert2edi(stationname.upper(),
                                     datadir,
                                     survey_configfile,
                                     birrp_configfile)

        except:
            raise
            print('Could not generate EDI file')
            edifn = None

    return edifn, cohfn

if __name__ == '__main__':
    main()
