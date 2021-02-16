# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby, YG


Plot RMS at each station as a map

"""
import os

from mtpy.modeling.modem import PlotRMSMaps
from tests import SAMPLE_DIR
from tests.imaging import ImageTestCase


class Test_PlotRMSMap(ImageTestCase):
    def test_fun(self):
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
        rmsmap = PlotRMSMaps(
            residual_fn=os.path.join(wd, filestem + ".res"),
            period_index=period_index,
            xminorticks=50000,
            yminorticks=50000,
            save_plots="y",
            plot_yn="n",
        )
        rmsmap.plot()

        rmsmap.save_figure(save_path, fig_close=False)  # this will save a file to
