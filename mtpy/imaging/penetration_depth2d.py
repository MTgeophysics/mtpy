"""
Description:
    With an input edi_file_folder and a list of period index,
    generate a profile using occam2d module,
    then plot the Penetration Depth profile at the given periods vs the stations locations.

Usage:
    python mtpy/imaging/penetration_depth2d.py /path2/edi_files_dir/   period_index_list
    python mtpy/imaging/penetration_depth2d.py.py examples/data/edi2/ 0 1 10 20 30 40

Author: fei.zhang@ga.gov.au
Date:   2017-01-23

Revision History:
    brenainn.moushall@ga.gov.au 03-04-2020 15:41:39 AEDT:
        - Modify 2D plot profile to take a list of selected periods
          instead of period indicies
"""

import glob
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import click

from mtpy.imaging.penetration import get_penetration_depth_by_index, load_edi_files, Depth2D

# mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'

# mpl.rcParams['figure.figsize'] = [20, 10]

import mtpy.core.mt as mt
import mtpy.modeling.occam2d_rewrite as occam2d_new
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
_logger = MtPyLog.get_mtpy_logger(__name__)
# logger =
# MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) #
# specific


# use the Zcompotent=[det, zxy, zyx]
def plot2Dprofile(edi_dir, selected_periods, ptol=0.05, zcomponent='det',
                  edi_list=None, tick_params={}, save=False, savepath=None, **kwargs):
    edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))

    _logger.debug("edi files: %s", edifiles)

    edis = load_edi_files(edi_dir, file_list=edi_list)
    plot = Depth2D(selected_periods, edis, ptol, zcomponent)
    plot.plot(tick_params, **kwargs)
    if save:
        if os.path.isdir(savepath):
            savepath == os.path.join(savepath, 'Depth2D.png')
        if savepath is not None:
            plot._fig.savefig(savepath)
        else:
            savepath = os.path.join(edi_dir, 'Depth2D.png')
    plot.show()


def barplot_multi_station_penentration_depth(
        edifiles_dir, per_index=0, zcomponent='det'):
    """
    A simple bar chart plot of the penetration depth across multiple edi files (stations),
    at the given (frequency) per_index. No profile-projection is done in this funciton.
    :param edifiles_dir: a list of edi files, or a dir of edi
    :param per_index: an integer smaller than the number of MT frequencies in the edi files.
    :return:
    """

    if os.path.isdir(edifiles_dir):
        edi_dir = edifiles_dir  # "E:/Githubz/mtpy2/tests/data/edifiles/"
        edifiles_dir = glob.glob(os.path.join(edi_dir, '*.edi'))
        _logger.debug(edifiles_dir)
    else:
        # Assume edifiles_dir is [a list of edi files]
        pass

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    # per_index=0,1,2,....
    periods = []

    depths = []

    stations = []

    mt_obj_list = [mt.MT(afile) for afile in edifiles_dir]

    (stations, periods, depths, _) = get_penetration_depth_by_index(
        mt_obj_list, int(per_index), whichrho=zcomponent)

    # the attribute Z
    # zeta = mt_obj.Z
    #
    # if per_index >= len(zeta.freq):
    #     raise Exception("Error: input period index must be less than number of freqs in zeta.freq=%s",len(zeta.freq))
    #
    # per = 1.0 / zeta.freq[per_index]
    # periods.append(per)
    # penetration_depth = -scale_param * np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
    #
    # depths.append(penetration_depth)
    # stations.append(mt_obj.station)

    #plt.plot(app_resis, color='b', marker='o')

    index = np.arange(len(depths))

    plt.bar(index, depths, color='#000000')

    # plt.xaxis.tick_top()
    # plt.set_xlabel('X LABEL')
    # plt.xaxis.set_label_position('top')

    plt.xlabel(
        'Penetration Depth Across Stations, for MT period= %6.5f Seconds' %
        periods[0], fontsize=16)
    plt.ylabel('Penetration Depth (m)', fontsize=16)
    # plt.title('Penetration Depth profile for T=??')
    bar_width = 0.4
    plt.xticks(
        index + bar_width / 2,
        stations,
        rotation='horizontal',
        fontsize=14)
    plt.legend()

    # plt.tight_layout()
    plt.gca().xaxis.tick_top()
    plt.show()

    # Check that the periods are the same value for all stations
    return (stations, depths, periods)


# =============================================================================================
# Example Usage:
# python mtpy/imaging/penetration_depth2d.py examples/data/edi_files/ 1 10 20 30
# python mtpy/imaging/penetration_depth2d.py tests/data/edifiles/ 0 1 10 20 30 40 50 59
# python mtpy/imaging/penetration_depth2d.py examples/data/edi2/ 0 1 10 20 30 40
# =============================================================================================
if __name__ == "__main__old":

    if len(sys.argv) < 2:
        print(("Usage: %s edi_dir" % sys.argv[0]))
        print("python examples/penetration_depth2d.py tests/data/edifiles/ 0 1 10 20 30 40 50 59")
        sys.exit(1)
    elif os.path.isdir(sys.argv[1]):
        edi_dir = sys.argv[1]  # the first argument is path2_edi_dir
        # the second,.... will be period index list
        period_index_list = sys.argv[2:]
        print(("period_index_list = {}".format(period_index_list)))

        # the rho zcomponent can be det, zxy zyx
        plot2Dprofile(edi_dir, period_index_list, zcomponent='det')
        perindex = int(sys.argv[2])
        # barplot_multi_station_penentration_depth(edi_dir, per_index=perindex)
        # #, zcomponent='zxy')
    else:
        print("Please provide an edi directory and period_index_list")


# =============================================================================================
# Command line wrapper
# =============================================================================================

@click.command()
@click.option('-i', '--input', type=str, default='examples/data/edi_files', help='directory or edsi data files')
@click.option('-p', '--period_list', type=str, default="0 1 10 20 30 40", help='Periods seperated by space')
def plot_penetration_image(input, period_list):
    if os.path.isdir(input):
        period_index_list = period_list.split(' ')
        plot2Dprofile(input, period_index_list, zcomponent='det')
    else:
        print("Please provide an edi directory !")


if __name__ == '__main__':
    plot_penetration_image()
