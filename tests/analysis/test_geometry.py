# -*- coding: utf-8 -*-
"""
TEST mtpy.core.mt.MT

@author: YG
"""
from unittest import TestCase
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:19:35 2017

@author: u64125
Fei.zhang@ga.gov.au
YG
"""

import os
import numpy as np

from mtpy.core.mt import MT
from tests import TEST_MTPY_ROOT
import mtpy.analysis.geometry as mtg


class Test_Geometry(TestCase):
    def test_get_dimensionality_from_edi_file(self):
        mt_obj = MT(os.path.normpath(os.path.join(TEST_MTPY_ROOT, "examples/data/edi_files/pb42c.edi")))
        dimensionality_result = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                                          2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        dimensionality = mtg.dimensionality(z_object=mt_obj.Z)

        self.assertTrue(np.allclose(dimensionality, dimensionality_result, 1e-8))

    def test_get_eccentricity_from_edi_file(self):
        mt_obj = MT(os.path.normpath(os.path.join(TEST_MTPY_ROOT, "examples/data/edi_files/pb42c.edi")))
        eccentricity_pb42c = (np.array([0.01675639, 0.01038589, 0.00527011, 0.00638819, 0.01483804,
                                        0.00385233, 0.00513294, 0.00403781, 0.02862114, 0.02689821,
                                        0.01425044, 0.05686524, 0.05742524, 0.02696736, 0.0275285,
                                        0.03647819, 0.04721932, 0.06336521, 0.12789841, 0.16409303,
                                        0.20630821, 0.34261225, 0.3967886, 0.51629705, 0.56645987,
                                        0.52558696, 0.46954261, 0.48028767, 0.47490701, 0.44927612,
                                        0.45185046, 0.44143159, 0.43570377, 0.41537978, 0.40546014,
                                        0.38785478, 0.37174031, 0.34534557, 0.35510941, 0.32282644,
                                        0.28501461, 0.22463964, 0.20683855]),
                              np.array([0.17132216, 0.2757994, 0.71263216, 0.50481657, 0.21604906,
                                        0.98931454, 0.75816349, 1.06885049, 0.23412284, 0.25015825,
                                        0.4117732, 0.06824775, 0.13024193, 0.49471091, 0.61126932,
                                        0.5471021, 0.6073574, 0.50578334, 0.30809787, 0.44938001,
                                        0.35430928, 0.20402482, 0.36750578, 0.30360427, 0.27660847,
                                        0.55139247, 0.53103062, 0.48771581, 0.19105325, 0.68542871,
                                        0.66189643, 0.1495947, 0.11353391, 0.09190586, 0.09006473,
                                        0.1079376, 0.13673274, 0.19349474, 0.23780856, 0.35159944,
                                        0.55386034, 0.78687532, 0.9654131]))

        eccentricity = mtg.eccentricity(z_object=mt_obj.Z)
        for i in range(2):
            self.assertTrue(np.allclose(eccentricity[i], eccentricity_pb42c[i], 1e-8))

    def test_get_strike_from_edi_file(self):
        edifile = os.path.normpath(os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files/pb42c.edi'))
        mt_obj = MT(edifile)
        strike_angle_pb42c = np.array([[np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [np.nan, np.nan],
                                       [38.45662316, 128.45662316],
                                       [28.61883115, 118.61883115],
                                       [14.45341494, 104.45341494],
                                       [8.43320651, 98.43320651],
                                       [4.94952784, 94.94952784],
                                       [2.09090369, 92.09090369],
                                       [1.39146887, 91.39146887],
                                       [0.39905337, 90.39905337],
                                       [-5.49553673, 84.50446327],
                                       [-6.28846049, 83.71153951],
                                       [-7.31641788, 82.68358212],
                                       [-10.45341947, 79.54658053],
                                       [-7.07075086, 82.92924914],
                                       [-7.5429295, 82.4570705],
                                       [-6.06405688, 83.93594312],
                                       [-3.54915951, 86.45084049],
                                       [-3.12596637, 86.87403363],
                                       [-0.47404093, 89.52595907],
                                       [2.74343665, 92.74343665],
                                       [4.78078759, 94.78078759],
                                       [7.71125988, 97.71125988],
                                       [11.0123521, 101.0123521],
                                       [13.81639678, 103.81639678],
                                       [13.60497071, 103.60497071],
                                       [15.87672806, 105.87672806]])

        strike_angle = mtg.strike_angle(z_object=mt_obj.Z)

        self.assertTrue(
            np.allclose(
                strike_angle[np.isfinite(strike_angle)],
                strike_angle_pb42c[np.isfinite(strike_angle_pb42c)],
                1e-8)
        )
