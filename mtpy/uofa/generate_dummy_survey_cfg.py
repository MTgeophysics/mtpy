#!/usr/bin/env python
"""
Convenience script for generating a dummy survey configuration file.

Useful for syntax testing or for synthetic data sets (e.g. converting BIRRP output to EDI).

arguments:
list of station names, separated by commas (no whitespaces!)

output:
dummy survey config file - 'dummy_survey.cfg'

"""

import os
import sys
import os.path as op
import mtpy.utils.filehandling as MTfh
import mtpy.utils.configfile as MTcf
#reload(MTfh)
#reload(MTcf)


def main():

    args = sys.argv[1:]

    if len(args) == 0:
        print('\n\tError - no station name given\n')
        return

    lo_stations = args[0]
    try:
        lo_stations = lo_stations.split(',')
        lo_stations = [i.upper() for i in lo_stations]
    except:
        print('Error - station name list could not be read!')
        print('(Must be comma-separated list...no whitespaces)\n')

    outfile = op.abspath(op.join(os.curdir, 'dummy_survey.cfg'))
    outfile = MTfh.make_unique_filename(outfile)

    survey_dict = {}

    for st in lo_stations:
        tmp_dict = {}

        for idx_r, req in enumerate(MTcf.list_of_required_keywords):
            tmp_dict[req] = MTcf.list_of_keyword_defaults_general[idx_r]
        for idx_e, efield in enumerate(MTcf.list_of_efield_keywords):
            tmp_dict[efield] = MTcf.list_of_keyword_defaults_efield[idx_e]
        for idx_b, bfield in enumerate(MTcf.list_of_bfield_keywords):
            tmp_dict[bfield] = MTcf.list_of_keyword_defaults_bfield[idx_b]

        survey_dict[st] = tmp_dict

    MTcf.write_dict_to_configfile(survey_dict, outfile)


if __name__ == '__main__':
    main()
