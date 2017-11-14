# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot responses from ModEM model.
"""
import os
import os.path as op
from unittest import TestCase

from mtpy.modeling.ModEM import PlotResponse

from tests import plt_wait, plt_close, SAMPLE_DIR
from tests.imaging import reset_matplotlib


class Test_ModEM_PlotResponse(TestCase):
    @classmethod
    def setUpClass(cls):
        reset_matplotlib()

    def tearDown(self):
        plt_wait(5)
        plt_close()

    def test_modular_MPI_NLCG_004(self):
        wd = op.normpath(op.join(SAMPLE_DIR, 'ModEM'))
        filestem = 'Modular_MPI_NLCG_004'
        datafn = 'ModEM_Data.dat'
        station = 'pb23'
        plot_z = False

        ro = PlotResponse(data_fn=op.join(wd, datafn),
                          resp_fn=op.join(wd, filestem + '.dat'),
                          plot_type=[station],
                          plot_z=plot_z)

        ro.plot()
