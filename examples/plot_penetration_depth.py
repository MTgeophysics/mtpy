"""
Description:
    Plot the Penetration Depth as a function of period (1/freq).
    The purpose is to show a profile over a collection of MT ground stations.

Author: fei.zhang@ga.gov.au

Date: 2017-01-23
"""

__author__ = 'fei.zhang@ga.gov.au'

import os
import sys
import glob
import matplotlib.pyplot as plt
import mtpy.imaging.plotresponse as mtpr
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
#logger = MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) # specific

def plot_edi_dir(edi_path):
    """ plot edi files from the input directory edi_path
    """

    edi_files = glob.glob(os.path.join(edi_path, "*.edi"))

    logger.debug(edi_files)

    for efile in edi_files:
    #for efile in edi_files[:2]:
        logger.debug("plotting %s", efile)
        # eo = mtedi.Edi(filename=efile)
        plot1(efile)

    return


def plot1(edi_file):
    """
    Plot the input edi_file
    Args:
        edi_file: path2edifile

    Returns:

    """

    # plt.style.use('dark_background')
    plt.style.use('seaborn-deep')
    plt.style.use('classic')

    logger.info("Plotting the edi file %s", edi_file)

    pr = mtpr.PlotResponse(fn=edi_file, plot_num=2, res_limits=(1, 10000), phase_limits=(0, 90))

    return pr


#########################################################
# plot one-by-one edi files in a given dirpath
# How to Run:
# export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH
# python plot_edis.py data/edi_files/

if __name__ == '__main__':
    edi_path = sys.argv[1]

    if os.path.isdir(edi_path):
        plot_edi_dir(edi_path)
    elif os.path.isfile(edi_path):
        plot1(edi_path)
    else:
        logger.error("Usage %s %s", sys.argv[0], "path2edi")
