# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:59:50 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

from mtpy import MT
from mtpy.core.mt_dataframe import MTDataFrame

from mt_metadata import TF_EDI_CGG

# =============================================================================


class TestMT(unittest.TestCase):
    def setUp(self):
        self.mt = MT()

    def test_clone_empty(self):
        self.mt.station = "test_01"
        self.mt.survey = "big"
        self.mt.latitude = 10
        self.mt.longitude = 20
        new_mt = self.mt.clone_empty()

        for attr in ["survey", "station", "latitude", "longitude"]:
            with self.subTest(attr):
                self.assertEqual(getattr(new_mt, attr), getattr(self.mt, attr))

        with self.subTest("tf is empty"):
            self.assertFalse(new_mt.has_transfer_function())

    def test_set_z(self):
        pass


class TestMT2DataFrame(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.m1 = MT(TF_EDI_CGG)
        self.m1.read_tf_file()

        self.sdf = self.m1.to_dataframe()
        self.mt_df = MTDataFrame()

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

    def test_to_z(self):
        self.assertEqual(self.m1.Z, self.mt_df.to_z_object(self.sdf))

    def test_to_t(self):
        self.assertEqual(self.m1.Tipper, self.mt_df.to_t_object(self.sdf))


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
