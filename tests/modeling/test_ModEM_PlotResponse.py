# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot responses from ModEM model.
"""
import os
import os.path as op
from unittest import TestCase

from mtpy.modeling.modem import PlotResponse

from tests.beta import SAMPLE_DIR
from tests import _plt_wait, _plt_close


class Test_ModEM_PlotResponse(TestCase):
    def tearDown(self):
        _plt_wait(5)
        _plt_close()

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
