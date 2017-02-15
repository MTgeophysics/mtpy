import sys
import os
import glob
import numpy as np


import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'

mpl.rcParams['figure.figsize']=[30,10]

import mtpy.core.mt as mt
import mtpy.modeling.occam2d_rewrite as occam2d_new
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
#logger = MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) # specific

def plot2Dprofile(edi_dir, period_index_list=None):
    #edi_dir = "/Softlab/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir="E:/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir=r"E:\Githubz\mtpy2\examples\data/edi2"

    #1 get a list of edi files, which are suppose to be in a profile.
    edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))

    logger.debug("edi files: ", edifiles)

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

    for period_index in period_index_list:   #[0,10,20,30]:

        logger.debug("doing period index %s", period_index)

        (stations, pen, periods)= get_penetration_depth(int(period_index), pr.edi_list)

        line_label="Period=%s s"%periods[0]

        plt.plot(pr.station_locations, pen, "--", marker='o', markersize="12", linewidth="2", label=line_label)
        plt.legend()

    plt.ylabel('Penetration Depth (m)', fontsize=16)
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


def plot_3D_profile(edi_dir, period_index):

    return "image"



def get_penetration_depth(per_index, mt_obj_list):

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

        penetration_depth = - scale_param * np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
        pen_depth.append(penetration_depth)

        stations.append(mt_obj.station)

    return(stations, pen_depth, periods)

# =============================================================================================
# python examples/plot_penetration_depth_profiles.py tests/data/edifiles/ 0 1 10 20 30 40 50 59
# python examples/plot_penetration_depth_profiles.py examples/data/edi2/ 0 1 10 20 30 40
# =============================================================================================
if __name__=="__main__":

    if len(sys.argv)<2:
        print("Usage: %s edi_dir"%sys.argv[0])
        print ("python examples/plot_penetration_depth_profiles.py tests/data/edifiles/ 0 1 10 20 30 40 50 59")

        sys.exit(1)
    elif os.path.isdir(sys.argv[1]):
        edi_dir = sys.argv[1]
        period_index_list=sys.argv[2:]
        plot2Dprofile(edi_dir, period_index_list )

    else:
        print("Please provide an edi directory and period_index_list")