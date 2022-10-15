# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:59:50 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

from mtpy.core.mt_dataframe import StationDataFrame
from mtpy import MT

from mt_metadata import TF_EDI_CGG

# =============================================================================


class TestStationDataFrame(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.sdf = StationDataFrame()
        m1 = MT(TF_EDI_CGG)
        m1.read_tf_file()

        self.sdf.from_tf(m1)

    def test_station(self):
        self.assertEqual(self.sdf.station, "TEST01")

    def test_period(self):
        self.assertEqual(self.sdf.period.size, 73)

    def test_frequency(self):
        self.assertEqual(self.sdf.frequency.size, 73)

    def test_latitude(self):
        self.assertEqual(self.sdf.latitude, -30.930285)

    def test_longitude(self):
        self.assertEqual(self.sdf.longitude, 127.22923)

    def test_elevation(self):
        self.assertEqual(self.sdf.elevation, 175.27)


class TestSetStationDataFrame(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.sdf = StationDataFrame()
        m1 = MT(TF_EDI_CGG)
        m1.read_tf_file()

        self.sdf.from_tf(m1)

    def test_set_station(self):
        self.sdf.station = "mt01"
        self.assertEqual(self.sdf.station, "mt01")

    def test_set_latitude(self):
        with self.subTest("string"):
            self.sdf.latitude = "40:00:00"
            self.assertEqual(self.sdf.latitude, 40.0)

        with self.subTest("float"):
            self.sdf.latitude = 50.0
            self.assertEqual(self.sdf.latitude, 50.0)

    def test_set_longitude(self):
        with self.subTest("string"):
            self.sdf.longitude = "-120:00:00"
            self.assertEqual(self.sdf.longitude, -120.0)

        with self.subTest("float"):
            self.sdf.longitude = -115.0
            self.assertEqual(self.sdf.longitude, -115.0)

    def test_set_elevation(self):
        self.sdf.elevation = 100
        self.assertEqual(self.sdf.elevation, 100)


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
