# -*- coding: utf-8 -*-
"""
plots phase tensor ellipses as a pseudo section (distance along profile vs period)

CreatedOn:      Wed Sep 18 15:35:39 2013
CreatedBy:      Alison Kirkby

LastUpdated:    2017-01-24
UpdatedBy:      fei.zhang@ga.gov.au

LastUpdated:    2017-11-24  FZ fixed this script after the big merge brokeness

"""


import os
import sys
import glob
from mtpy.core import mt
import click

import matplotlib.pyplot as plt
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
_logger = MtPyLog.get_mtpy_logger(__name__)
# logger =
# MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) #
# specific


def plot_edi_dir(edi_path, every_how_many_edi=2):
    """ plot edi files from the input directory edi_path
    """

    edi_files = glob.glob(os.path.join(edi_path, "*.edi"))

    _logger.debug(edi_files)

    for efile in edi_files[::every_how_many_edi]:
        # for efile in edi_files[:2]:
        _logger.debug("plotting %s", efile)
        # eo = mtedi.Edi(filename=efile)
        plot_edi_file(efile)

        # plt.pause(1.0)
        # plt.close()

    return


def plot_edi_file(edi_file):
    """
    Plot the input edi_file
    Args:
        edi_file: path2edifile

    Returns:

    """

    # plt.style.use('dark_background')
    plt.style.use('seaborn-deep')
    plt.style.use('classic')

    _logger.info("Plotting the edi file %s", edi_file)

    mt_obj = mt.MT(edi_file)
    pt_obj = mt_obj.plot_mt_response(plot_yn='n')
    pt_obj.plot()

    # pt_obj = mt_obj.plot_mt_response(plot_yn='n',plot_num=2, res_limits=(1, 10000), phase_limits=(0, 90))
    # pt_obj.plot()

    return



@click.command()
@click.option('-d','--directory',type=str,default='examples\data\edi_files',help='directory of edsi files (file/directory')
@click.option('-c','--count',type=int,default=6, help='every how many edsi files')
def select_plot_edi_files(directory,count):
    print('Directory of edsi files ------> {}'.format(directory))
    print('Count of files          ------> {}'.format(count))
    if os.path.isfile(directory):
        plot_edi_file(directory)
    elif os.path.isdir(directory):
        #plot_edi_dir(edi_path)
        # plot_edi_dir(edi_path,every_how_many_edi=6)
        plot_edi_dir(directory,every_how_many_edi=count)

###############################################################################
# plot one-by-one edi files in a given dir-path
# How to Run:
#    export PYTHONPATH=.
#    python  examples/plot_edis.py data/edifiles/
#    python  examples/plot_edis.py data/edifiles/15125A.edi
# =============================================================================
if __name__ == '__main__old':

    select_plot_edi_files()

    if len(sys.argv) < 2:
        print (
            "\n please provide path to edi files\n USAGE:  %s path2edifile" %
            sys.argv[0])
        sys.exit(1)
    else:
        edi_path = sys.argv[1]

        if os.path.isfile(edi_path):
            plot_edi_file(edi_path)
        elif os.path.isdir(edi_path):
            #plot_edi_dir(edi_path)
            plot_edi_dir(edi_path,every_how_many_edi=6)
        else:
            _logger.error("Usage %s %s", sys.argv[0], "path2edi")
###############################################################################
# Command wrapper for the edsi files plotting
# How to Run:
#           python  examples/scripts/plot_edis.py examples/data/edifiles
#           python  examples/scripts/plot_edis.py --help
###############################################################################

if __name__ == '__main__':

    select_plot_edi_files()


