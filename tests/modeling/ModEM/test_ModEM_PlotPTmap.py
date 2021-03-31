# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot data and responses from ModEM model.

"""

import os.path as op

# from mtpy.modeling.modem import PlotPTMaps
from mtpy.modeling.modem.phase_tensor_maps import PlotPTMaps
from tests import SAMPLE_DIR
from tests.imaging import ImageTestCase


class Test_ModEM_PlotPTMaps(ImageTestCase):
    def test_modular_MPI_NLCG_004(self):
        wd = op.normpath(op.join(SAMPLE_DIR, "ModEM"))
        filestem = "Modular_MPI_NLCG_004"
        datafn = "ModEM_Data.dat"
        PlotPTMaps(
            data_fn=op.join(wd, datafn),
            resp_fn=op.join(wd, filestem + ".dat"),
            ellipse_size=20,
        )
