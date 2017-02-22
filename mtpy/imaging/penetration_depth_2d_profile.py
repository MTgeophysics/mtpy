"""
Description:
    With an input edi_file_folder and a list of period index,
    generate a profile using occam2d module,
    then plot the Penetration Depth profile at the given periods vs the stations locations.

Usage:
    python mtpy/imaging/penetration_depth_2d_profile.py /path2/edi_files_dir/   period_index_list
    python mtpy/imaging/penetration_depth_2d_profile.py examples/data/edi2/ 0 1 10 20 30 40

Author: fei.zhang@ga.gov.au
Date:   2017-01-23
"""

import sys
import os
import glob
import numpy as np


import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'

mpl.rcParams['figure.figsize']=[20,10]

import mtpy.core.mt as mt
import mtpy.modeling.occam2d_rewrite as occam2d_new
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
#logger = MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) # specific

def plot2Dprofile(edi_dir, period_index_list=None, zcomponent='det'): #use the Zcompotent=[det, zxy, zyx]
    #edi_dir = "/Softlab/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir="E:/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir=r"E:\Githubz\mtpy2\examples\data/edi2"

    #1 get a list of edi files, which are suppose to be in a profile.
    edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))

    logger.debug("edi files: %s", edifiles)

    # stations = ['151{0:02}A'.format(s) for s in range(24, 31)]
    # pr = occam2d_new.Profile(edi_path=edi_dir, station_list=stations)

    pr = occam2d_new.Profile(edi_path=edi_dir)

    pr.generate_profile()

    pr.plot_profile(station_id=[0, 4])


    #2 get the period_index list

    if period_index_list == None or len(period_index_list)==0:
        logger.error("Please provide a period index list like [1,2,3,4]")
        raise Exception("Period index list is empty")

    #  set station labels to only be from 1st to 4th index of station name
    plt.figure()
    for period_index in period_index_list:   #[0,10,20,30]:

        logger.debug("doing period index %s", period_index)

        (stations, pen, periods)= get_penetration_depth(pr.edi_list, int(period_index), whichrho=zcomponent)

        line_label="Period=%s"%periods[0]
        
        plt.plot(pr.station_locations, pen, "--", marker='o', markersize="12", linewidth="2", label=line_label)
        plt.legend()

    plt.ylabel('Penetration Depth (Metres) Computed by %s'%zcomponent, fontsize=16)
    plt.yticks(fontsize=16)

    plt.xlabel('MT Penetration Depth Profile Over Stations.', fontsize=16)
    logger.debug("stations= %s", stations)
    logger.debug("station locations: %s",pr.station_locations)
    if (pr.station_list is not None):
        plt.xticks( pr.station_locations, pr.station_list, rotation='horizontal', fontsize=16)
    else:  # Are the stations in the same order as the profile generated pr.station_list????
        plt.xticks(pr.station_locations, stations, rotation='horizontal', fontsize=16)

    # plt.tight_layout()
    plt.gca().xaxis.tick_top()

    plt.show()


def get_penetration_depth(mt_obj_list, per_index, whichrho='det'): #whichrho=[det, zxy, zyx]
    """
    compute the penetration depth of mt_obj at the given period_index, and using whichrho option
    :param per_index:
    :param mt_obj_list:
    :param whichrho:
    :return:
    """

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    # per_index=0,1,2,....
    periods = []

    pen_depth= []

    stations = []

    for mt_obj in mt_obj_list:

        # the attribute Z
        zeta = mt_obj.Z

        if per_index >= len(zeta.freq):
            logger.debug ("number of frequecies= %s", len(zeta.freq) )
            raise Exception("Index out_of_range Error: period index must be less than number of periods in zeta.freq")

        per = 1.0 / zeta.freq[per_index]
        periods.append(per)

        if whichrho=='zxy':
            penetration_depth = - scale_param * np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
        elif whichrho=='zyx':
            penetration_depth = - scale_param * np.sqrt(zeta.resistivity[per_index, 1, 0] * per)
        elif whichrho =='det':
            # determinant
            det2 = np.abs(zeta.det[0][per_index])  # determinant value at the given period index
            penetration_depth = -scale_param * np.sqrt(0.2 * per * det2 * per)
        else:
            logger.critical("unsupported method to compute penetratoin depth: %s", whichrho)
            sys.exit(100)

        pen_depth.append(penetration_depth)

        stations.append(mt_obj.station)

    return(stations, pen_depth, periods)


def barplot_multi_station_penentration_depth(edifiles_dir, per_index=0, zcomponent='det'):
    """
    A simple bar chart plot of the penetration depth across multiple edi files (stations),
    at the given (frequency) per_index. No profile-projection is done in this funciton.
    :param edifiles_dir: a list of edi files, or a dir of edi
    :param per_index: an integer smaller than the number of MT frequencies in the edi files.
    :return:
    """

    if os.path.isdir(edifiles_dir):
        edi_dir = edifiles_dir # "E:/Githubz/mtpy2/tests/data/edifiles/"
        edifiles_dir = glob.glob(os.path.join(edi_dir, '*.edi'))
        logger.debug(edifiles_dir)
    else:
        # Assume edifiles_dir is [a list of edi files]
        pass

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    # per_index=0,1,2,....
    periods = []

    depths = []

    stations = []

    mt_obj_list=[mt.MT(afile) for afile in edifiles_dir]

    (stations, depths, periods) = get_penetration_depth(mt_obj_list, int(per_index), whichrho=zcomponent)

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

    plt.xlabel('Penetration Depth Across Stations, for MT period= %6.5f Seconds' % periods[0], fontsize=16)
    plt.ylabel('Penetration Depth (m)', fontsize=16)
    # plt.title('Penetration Depth profile for T=??')
    bar_width = 0.4
    plt.xticks(index + bar_width / 2, stations, rotation='horizontal', fontsize=14)
    plt.legend()

    # plt.tight_layout()
    plt.gca().xaxis.tick_top()
    plt.show()

    #Check that the periods are the same value for all stations
    return (stations, depths, periods)

# =============================================================================================
# Example Usage:
# python mtpy/imaging/penetration_depth_2d_profile.py examples/data/edi_files/ 1 10 20 30
# python mtpy/imaging/penetration_depth_2d_profile.py tests/data/edifiles/ 0 1 10 20 30 40 50 59
# python mtpy/imaging/penetration_depth_2d_profile.py examples/data/edi2/ 0 1 10 20 30 40
# =============================================================================================
if __name__=="__main__":

    if len(sys.argv)<2:
        print("Usage: %s edi_dir"%sys.argv[0])
        print ("python examples/penetration_depth_2d_profile.py tests/data/edifiles/ 0 1 10 20 30 40 50 59")
        sys.exit(1)
    elif os.path.isdir(sys.argv[1]):
        edi_dir = sys.argv[1] # the first argument is path2_edi_dir
        period_index_list=sys.argv[2:] # the second,.... will be period index list

        plot2Dprofile(edi_dir, period_index_list, zcomponent='det') # the rho zcomponent can be det, zxy zyx
        perindex = int(sys.argv[2])
        #barplot_multi_station_penentration_depth(edi_dir, per_index=perindex) #, zcomponent='zxy')
    else:
        print("Please provide an edi directory and period_index_list")
