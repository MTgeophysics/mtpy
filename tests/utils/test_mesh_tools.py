# -*- coding: utf-8 -*-
"""

Created on Mon Nov 29 17:20:47 2021

:author: Jared Peacock

:license: MIT

"""

import unittest
import numpy as np

from mtpy.utils import mesh_tools

class TestMeshUtils(unittest.TestCase):
    
    def setUp(self):
        
        self.cell_sizes = [.1, .3, 1.0, 3, 10, 30, 100, 300, 1000, 3000]
        
    def test_rounding(self):
        for width in self.cell_sizes:
            with self.subTest(msg=f"width {width}"):
                rounding = int(-1 * np.floor(np.log10(width)))
                self.assertEquals(rounding, mesh_tools.get_rounding(width))
                
    def test_padding(self):
        pad = np.array([  1.,   3.,   4.,   6.,   9.,  14.,  23.,  37.,  61., 100.])
        self.assertTrue((pad == mesh_tools.get_padding_cells(1, 100, 10, 1.2)).all())
        
    def test_stretch(self):
        pad = np.array([ 1.,  2.,  3.,  5.,  7.,  9., 12., 16., 20., 25.])
        self.assertTrue((pad == mesh_tools.get_padding_from_stretch(1, 1.2, 10)).all())
                
    
                
if "main" in __name__:
    unittest.main()               