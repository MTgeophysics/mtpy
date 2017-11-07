# -*- coding: utf-8 -*-
"""
TEST mtpy.core.mt.MT

@author: YG
"""
from unittest import TestCase

import os
import numpy as np

from mtpy.core.mt import MT
from tests import TEST_MTPY_ROOT
import mtpy.analysis.geometry as mtg


class Test_MT(TestCase):
    def test_get_dimensionality_from_edi_file(self):
        mt_obj = MT(os.path.normpath(os.path.join(TEST_MTPY_ROOT, "data/edi_files/pb42c.edi")))
        dimensionality_result = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                                          2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        dimensionality = mtg.dimensionality(z_object=mt_obj.Z)

        self.assertTrue(np.allclose(dimensionality, dimensionality_result))
