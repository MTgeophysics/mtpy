# -*- coding: utf-8 -*-
"""
Plot Depth Slice
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby


Fei.zhang@ga.gov.au
"""

import os
import shutil
import matplotlib.pyplot as plt

from mtpy.modeling.modem import PlotDepthSlice

from tests.beta import *

def test_fun():
    """
    test function
    :return: T/F
    """

    # directory where files are located
    wd = os.path.join(SAMPLE_DIR,'ModEM')  # r'E:\Githubz\mtpy\examples\model_files\ModEM'

    # directory to save to
    save_path =  TEMP_OUT_DIR  # r'E:\Githubz\mtpy\temp'

    # file stem for inversion result
    filestem = 'Modular_MPI_NLCG_004'

    # period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
    period_index = 0

    # plot map
    dsmap = PlotDepthSlice(model_fn = os.path.join(wd,filestem+'.rho'),
                            data_fn = os.path.join(wd,filestem+'dat'),
                            depth_index=30,
                            save_plots='n'
                            )

    path2file = os.path.join(save_path,'DepthSlice.png')

    if os.path.exists(path2file):
        os.remove(path2file)

    plt.savefig(path2file)

    assert (os.path.exists(path2file))

#########################################
if __name__ == "__main__":
    test_fun()
