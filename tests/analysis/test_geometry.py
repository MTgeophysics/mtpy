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
        dimensionality_result = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                                          2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3])
        dimensionality = mtg.dimensionality(z_object=mt_obj.Z)

        self.assertTrue(np.allclose(dimensionality, dimensionality_result, 1e-8))

    def test_get_eccentricity_from_edi_file(self):
        mt_obj = MT(os.path.normpath(os.path.join(TEST_MTPY_ROOT, "examples/data/edi_files/pb42c.edi")))
        eccentricity_pb42c = (np.array(
            [0.01675639, 0.01038589, 0.00527011, 0.00638819, 0.01483804,
             0.00385233, 0.00513294, 0.00403781, 0.02862114, 0.02689821,
             0.01425044, 0.05686524, 0.05742524, 0.02696736, 0.0275285 ,
             0.03647819, 0.04721932, 0.06336521, 0.12789841, 0.16409303,
             0.20630821, 0.34261225, 0.3967886 , 0.51629705, 0.56645987,
             0.52558696, 0.46954261, 0.48028767, 0.47490701, 0.44927612,
             0.45185046, 0.44143159, 0.43570377, 0.41537978, 0.40546014,
             0.38785478, 0.37174031, 0.34534557, 0.35510941, 0.32282644,
             0.28501461, 0.22463964, 0.20683855]),
            np.array(
                [1.09157801, 1.88513021, 2.52461839, 2.82831998, 2.48686079,
                 2.78766939, 2.37250674, 1.16767867, 2.82659526, 1.41909434,
                 1.57463881, 2.80007766, 2.73439937, 2.71321983, 2.82682667,
                 2.25051091, 1.435991  , 0.43231073, 0.81874286, 1.62593021,
                 2.4864899 , 2.91852175, 3.11049976, 3.48724331, 3.6630486 ,
                 3.54297703, 3.40317309, 3.4175019 , 3.3910377 , 3.28570849,
                 3.33766706, 3.31489955, 3.32710981, 3.29740796, 3.28179559,
                 3.24496124, 3.2052965 , 3.13248919, 3.11445279, 2.98672687,
                 2.85242595, 2.79044042, 2.73406595]))

        eccentricity = mtg.eccentricity(z_object=mt_obj.Z)
        for i in range(2):
            self.assertTrue(np.allclose(eccentricity[i], eccentricity_pb42c[i], 1e-8))

    def test_get_strike_from_edi_file(self):
        edifile = os.path.normpath(os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files/pb42c.edi'))
        mt_obj = MT(edifile)
        strike_angle_pb42c = np.array([[      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [ 38.45662316, -51.54337684],
                                       [ 28.61883115, -61.38116885],
                                       [ 14.45341494, -75.54658506],
                                       [  8.43320651, -81.56679349],
                                       [  4.94952784, -85.05047216],
                                       [  2.09090369, -87.90909631],
                                       [  1.39146887, -88.60853113],
                                       [  0.39905337, -89.60094663],
                                       [ -5.49553673,  84.50446327],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [      np.nan,       np.nan],
                                       [ -6.06405688,  83.93594312],
                                       [ -3.54915951,  86.45084049],
                                       [ -3.12596637,  86.87403363],
                                       [ -0.47404093,  89.52595907],
                                       [  2.74343665, -87.25656335],
                                       [  4.78078759, -85.21921241],
                                       [  7.71125988, -82.28874012],
                                       [ 11.0123521 , -78.9876479 ],
                                       [ 13.81639678, -76.18360322],
                                       [ 13.60497071, -76.39502929],
                                       [      np.nan,       np.nan]])


        strike_angle = mtg.strike_angle(z_object=mt_obj.Z)

        self.assertTrue(
            np.allclose(
                strike_angle[np.isfinite(strike_angle)],
                strike_angle_pb42c[np.isfinite(strike_angle_pb42c)],
                1e-8)
        )
