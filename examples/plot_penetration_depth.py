"""
Description:
    Plot the Penetration Depth as a function of period (1/freq).
    The purpose is to show a profile over a collection of MT ground stations.

Author: fei.zhang@ga.gov.au
Date:   2017-01-23
"""

import glob
import os
import sys
import numpy as np

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.color'] = 'r'

import matplotlib.pyplot as plt

import mtpy.core.mt as mt

import mtpy.imaging.plot_mt_response as mtpr
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to config the logger
logger = MtPyLog().get_mtpy_logger(__name__)

# logger = MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) # specific

def plot_edi_dir(edi_path,rholist=['zxy','zyx','det']):
    """ plot edi files from the input directory edi_path
    """

    edi_files = glob.glob(os.path.join(edi_path, "*.edi"))

    logger.debug(edi_files)

    for efile in edi_files:
        # for efile in edi_files[:2]:
        logger.debug("plotting %s", efile)
        # eo = mtedi.Edi(filename=efile)
        plot_edi_file(efile,rholist=rholist)

    return


def plot_edi_file(edifile, rholist=['zxy','zyx','det'], savefile=None):
    """
    Plot the input edi_file
    Args:
        edi_file: path2edifile
        rholist: a list of the rho to be used.
        savefile: path2savefig, not save if None
    Returns:
    """
    # plt.style.use('dark_background')
    # plt.style.use('seaborn-deep')
    # plt.style.use('classic')
    plt.grid(True)

    logger.info("Plotting the edi file %s", edifile)

    mt_obj = mt.MT(edifile)
    zeta = mt_obj.Z     # the attribute Z represent the impedance tensor 2X2 matrix
    freqs = zeta.freq   # frequencies

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    logger.debug("scale parameter= %s",scale_param)

    # The periods array

    periods = 1.0 / freqs

    legendh=[]

    if 'zxy' in rholist:
        # One of the 4-components: XY
        penetration_depth = scale_param * np.sqrt(zeta.resistivity[:, 0, 1] * periods)

        #pen_zxy, = plt.semilogx(periods, -penetration_depth, '-*',label='Zxy')
        pen_zxy, = plt.semilogx(periods, -penetration_depth, color='#000000', marker='*', label='Zxy')
        #See http://matplotlib.org/1.3.1/examples/pylab_examples/line_styles.html

        legendh.append(pen_zxy)

        # plt.title("XY Penetration Depth in Meters")
        # plt.xlabel("Period (seconds)")
        # plt.ylabel("Depth Meters")

    if 'zyx' in rholist:
        penetration_depth = scale_param * np.sqrt(zeta.resistivity[:, 1, 0] * periods)

        pen_zyx, = plt.semilogx(periods, -penetration_depth, color='g', marker='o', label='Zyx')
        legendh.append(pen_zyx)

    if 'det' in rholist:
        # determinant
        det2 = np.abs(zeta.det[0])
        det_penetration_depth = scale_param * np.sqrt(0.2 * periods * det2 * periods)

        #pen_det, = plt.semilogx(periods, -det_penetration_depth, '-^', label='Determinant')
        pen_det, = plt.semilogx(periods, -det_penetration_depth, color='b', marker='^', label='Determinant')
        legendh.append(pen_det)


    plt.legend(handles=legendh, bbox_to_anchor=(0.1,0.5), loc=3, ncol=1, borderaxespad=0.)

    plt.title("Penetration Depth for file %s"% edifile)

    plt.xlabel("Log Period (seconds)",fontsize=16)

    plt.ylabel("Penetration Depth (meters)",fontsize=16)

    if savefile is not None:
        plt.savefig(savefile, dpi=800)

    plt.show()

    return savefile


def plot_multi_station_pen_depth(per_index, edifiles):
    """
    plot the pen depth of multiple edi files (stations), at the given frequency defined by per_index=0,1,2...
    :param edifiles:
    :param per_index:
    :return:
    """

    if os.path.isdir(edifiles):
        edi_dir = edifiles # "E:/Githubz/mtpy2/tests/data/edifiles/"
        edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))
        logger.debug(edifiles)
    else:
        # Assume edifiles is [a list of files]
        pass

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    # per_index=0,1,2,....
    periods = []

    app_resis = []

    stations = []

    for afile in edifiles:
        mt_obj = mt.MT(afile)

        # the attribute Z
        zeta = mt_obj.Z

        if per_index >= len(zeta.freq):
            raise Exception("Error: input period index must be less than number of freqs in zeta.freq=%s",len(zeta.freq))

        per = 1.0 / zeta.freq[per_index]
        periods.append(per)
        penetration_depth = -scale_param * np.sqrt(zeta.resistivity[per_index, 0, 1] * per)

        app_resis.append(penetration_depth)
        stations.append(mt_obj.station)

    #plt.plot(app_resis, color='b', marker='o')

    index = np.arange(len(app_resis))

    plt.bar(index, app_resis, color='#000000')

    # plt.xaxis.tick_top()
    # plt.set_xlabel('X LABEL')
    # plt.xaxis.set_label_position('top')

    plt.xlabel('Penetration Depth Across Stations, for MT period= %6.3f Seconds' % periods[0], fontsize=16)
    plt.ylabel('Penetration Depth (m)', fontsize=16)
    # plt.title('Penetration Depth profile for T=??')
    bar_width = 0.4
    plt.xticks(index + bar_width / 2, stations, rotation='horizontal', fontsize=16)
    plt.legend()

    # plt.tight_layout()
    plt.gca().xaxis.tick_top()
    plt.show()

    #Check that the periods are the same value for all stations
    return (stations, app_resis, periods)

###############################################################################
# plot one-by-one edi files in a given dirpath
# How to Run:
#    export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH
#    python  examples/plot_penetration_depth.py data/edi_files/
#    python  examples/plot_penetration_depth.py tests/data/edifiles/15125A.edi

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print ("\n please provide path to edi files\n USAGE:  %s path2edifile" % sys.argv[0])
        sys.exit(1)
    else:
        edi_path = sys.argv[1]

        if os.path.isfile(edi_path):
            plot_edi_file(edi_path , savefile='C:/temp/pen_depth.jpg')
            # rholist can be any of ['zxy','zyx','det'], default all of them
        elif os.path.isdir(edi_path):
            #plot_edi_dir(edi_path )
            #plot_edi_dir(edi_path, rholist=['det'] )
            plot_multi_station_pen_depth(42, edi_path)
        else:
            logger.error("Usage %s %s", sys.argv[0], "path2edi")
