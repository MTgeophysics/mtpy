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

    def test_set_station(self):
        self.sdf.station = "mt01"
        self.assertEqual(self.sdf.station, "mt01")


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
