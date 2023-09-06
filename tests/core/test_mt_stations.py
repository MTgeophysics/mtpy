# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 16:27:01 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import unittest

from mtpy.core import MTStations
from mtpy import MT

# =============================================================================


class TestMTStationGrid(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.east = 243900.352
        self.north = 4432069.056898517
        self.utm_epsg = 32611
        dx = 1000
        dy = 2000
        count = 1
        mt_list = []
        for ii in range(5):
            for jj in range(5):
                mt_obj = MT(
                    east=(self.east + ii * dx),
                    north=(self.north + jj * dy),
                    utm_epsg=self.utm_epsg,
                    station=f"mt{count:02}",
                )
                mt_list.append(mt_obj)

        self.stations = MTStation(self.utm_epsg, mt_list=mt_list)

    def test_station_len(self):
        self.assertEqual(25, self.stations.station_locations.shape[0])


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
