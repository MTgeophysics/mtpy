# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot data and responses from ModEM model.

"""

import os.path as op
from unittest import TestCase

from mtpy.modeling.modem import PlotPTMaps
from tests.beta import SAMPLE_DIR

from tests import _plt_wait, _plt_close


class Test_ModEM_PlotPTMaps(TestCase):
    def tearDown(self):
        _plt_wait(5)
        _plt_close()

    def test_modular_MPI_NLCG_004(self):
        wd = op.normpath(op.join(SAMPLE_DIR, 'ModEM'))
        filestem = 'Modular_MPI_NLCG_004'
        datafn = 'ModEM_Data.dat'
        PlotPTMaps(data_fn=op.join(wd, datafn),
                   resp_fn=op.join(wd, filestem + '.dat'),
                   ellipse_size=20
                   )
