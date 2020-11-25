# -*- coding: utf-8 -*-
"""
Plot Depth Slice
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby


Fei.zhang@ga.gov.au
YG
"""

import os

import matplotlib.pyplot as plt

# from mtpy.modeling.modem import PlotDepthSlice
from mtpy.imaging.plot_depth_slice import PlotDepthSlice
from tests import SAMPLE_DIR
from tests.imaging import ImageTestCase


class Test_PlotDepthSlice(ImageTestCase):
    def test_PlotDepthSlice(self):
        """
        test function
        :return: T/F
        """

        # directory where files are located
        wd = os.path.join(SAMPLE_DIR, "ModEM")

        # directory to save to
        save_path = self._temp_dir
        # file stem for inversion result
        filestem = "Modular_MPI_NLCG_004"

        # period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
        period_index = 0

        # plot map
        dsmap = PlotDepthSlice(
            model_fn=os.path.join(wd, filestem + ".rho"),
            data_fn=os.path.join(wd, filestem + "dat"),
            depth_index=30,
            save_plots="n",
        )

        path2file = os.path.join(save_path, "DepthSlice.png")

        if os.path.exists(path2file):
            os.remove(path2file)

        plt.savefig(path2file)

        assert os.path.exists(path2file)
