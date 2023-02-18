# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:59:50 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

from mtpy.core.mt_dataframe import MTDataFrame
from mtpy import MT

from mt_metadata import TF_EDI_CGG

# =============================================================================


class TestMTDataFrame(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.m1 = MT(TF_EDI_CGG)
        self.m1.read_tf_file()

        self.sdf = self.m1.to_dataframe()

    def test_station(self):
        self.assertEqual(self.sdf.station.unique()[0], "TEST01")

    def test_period(self):
        self.assertEqual(self.sdf.period.size, 73)

    def test_latitude(self):
        self.assertEqual(self.sdf.latitude.unique()[0], -30.930285)

    def test_longitude(self):
        self.assertEqual(self.sdf.longitude.unique()[0], 127.22923)

    def test_elevation(self):
        self.assertEqual(self.sdf.elevation.unique()[0], 175.27)

    def test_to_z_object(self):
        new_z = self.sdf.to_z_object("TEST01")
        self.assertTrue(self.m1.Z == new_z)


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
