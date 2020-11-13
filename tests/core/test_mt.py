# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:32:44 2020

:author: Jared Peacock

:license: MIT

"""

import unittest
from mtpy.core import mt

class TestMT(unittest.TestCase):
    """
    test MT object functionality
    """
    
    def setUp(self):
        self.mt_obj = mt.MT()
        
    def test_position_ll(self):
        self.mt_obj.latitude = 40.0
        self.mt_obj.longitude = '-120:00:00'
        
        self.assertEqual(self.mt_obj.longitude, -120.0)
        self.assertEqual(self.mt_obj.latitude, 40.0)
        
        self.assertAlmostEqual(self.mt_obj.east, 243900, 0)
        self.assertAlmostEqual(self.mt_obj.north, 4432069, 0)
        self.assertEqual(self.mt_obj.utm_zone, "11S")
        
    def test_position_utm(self):
        self.mt_obj.east = 243900.35
        self.mt_obj.north = "4432069.0569"
        self.mt_obj.utm_zone = '11S'
        
        self.assertEqual(self.mt_obj.north, 4432069.0569)
        
        self.assertAlmostEqual(self.mt_obj.latitude, 40, 0)
        self.assertAlmostEqual(self.mt_obj.longitude, -120, 0)
        
        
        
        
# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
    