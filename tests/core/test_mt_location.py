# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:17:37 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

from mtpy.core.mt_location import MTLocation

# =============================================================================


class TestMTLocation(unittest.TestCase):
    def setUp(self):
        self.loc = MTLocation()

    def test_set_latitude(self):
        with self.subTest("string"):
            self.loc.latitude = "40:00:00"
            self.assertEquals(self.loc.latitude, 40.0)

        with self.subTest("int"):
            self.loc.latitude = 40
            self.assertIsInstance(self.loc.latitude, float)

        def set_lat(value):
            self.loc.latitude = value

        with self.subTest("Fails"):
            self.assertRaises(ValueError, set_lat, "ten")

        with self.subTest("Too large"):
            self.assertRaises(ValueError, set_lat, "100")

    def test_set_utm_zone(self):
        self.loc.latitude = 40
        self.loc.longitude = -120


# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    unittest.main()
