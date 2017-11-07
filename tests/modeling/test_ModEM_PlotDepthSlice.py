# -*- coding: utf-8 -*-
"""
Plot Depth Slice
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby


Fei.zhang@ga.gov.au
YG
"""
from unittest import TestCase

import os

import matplotlib.pyplot as plt
from mtpy.modeling.modem import PlotDepthSlice
from tests import TEST_TEMP_DIR
from tests.beta import SAMPLE_DIR
from tests.imaging import _plt_wait, _plt_close


class Test_PlotDepthSlice(TestCase):
    def tearDown(self):
        _plt_wait(1)
        _plt_close()

    def test_PlotDepthSlice(self):
        """
        test function
        :return: T/F
        """

        # directory where files are located
        wd = os.path.join(SAMPLE_DIR, 'ModEM')

        # directory to save to
        save_path = TEST_TEMP_DIR
        # file stem for inversion result
        filestem = 'Modular_MPI_NLCG_004'

        # period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
        period_index = 0

        # plot map
        dsmap = PlotDepthSlice(model_fn=os.path.join(wd, filestem + '.rho'),
                               data_fn=os.path.join(wd, filestem + 'dat'),
                               depth_index=30,
                               save_plots='n'
                               )

        path2file = os.path.join(save_path, 'DepthSlice.png')

        if os.path.exists(path2file):
            os.remove(path2file)

        plt.savefig(path2file)

        assert (os.path.exists(path2file))
