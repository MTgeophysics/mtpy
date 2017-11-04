# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot RMS at each station as a map

"""
import os

# not needed os.chdir(r'C:\Git\mtpy')
from mtpy.modeling.modem import Plot_RMS_Maps
from tests.beta import *

import matplotlib.pyplot as plt
plt.ion() # make figure disappear automatically

#plt.ioff()  # to make figure show normally, and require click to close the figure to continue

def test_fun():
    """
    test function
    :return: T/F
    """

    # directory where files are located
    wd =  os.path.join(SAMPLE_DIR,'ModEM')

    # directory to save to
    save_path = TEMP_OUT_DIR

    # file stem for inversion result
    filestem = 'Modular_MPI_NLCG_004'

    # period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
    period_index = 0

    # plot map
    rmsmap = Plot_RMS_Maps(residual_fn=os.path.join(wd,filestem + '.res'),period_index=period_index,
                       xminorticks=50000,yminorticks=50000,save_plots='y')

    rmsmap.save_figure(save_path) # this will save a file to E:\Githubz\mtpy\temp\00_RMS_0.020535_s..png

############################################
if __name__ == "__main__":
    test_fun()
