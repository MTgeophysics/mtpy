# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:44:40 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""

import unittest

import numpy as np

from mtpy.modeling.modem import Data
from tests import EDI_DATA_DIR, EDI_DATA_DIR2


class TestModEMData(unittest.TestCase):
    """
    test making a modem file
    """

    def setUp(self):
        self.data = Data(
            edi_list=list(EDI_DATA_DIR.glob("*.edi")),
            period_list=np.logspace(-3, -1, 20),
        )
        self.data.data_array = self.data.fill_data_array(self.data.mt_dict)

    def test_station_list(self):
        station_list = [fn.stem.replace("c", "") for fn in EDI_DATA_DIR.glob("*.edi")]

        self.assertListEqual(station_list, self.data.station_locations.station.tolist())

    def test_period_range(self):
        self.assertEqual(20, self.data.period_list.size)

    def test_add_station(self):
        edi_list = list(EDI_DATA_DIR2.glob("*.edi"))[0:2]
        new_data, new_mt_dict = self.data.add_station(fn=edi_list)

        self.assertIn(edi_list[0].stem, new_mt_dict.keys())
        self.assertIn(edi_list[1].stem, new_mt_dict.keys())
        self.assertIn(edi_list[0].stem, new_data["station"])
        self.assertIn(edi_list[1].stem, new_data["station"])

    def test_remove_station(self):
        new_data, new_mt_dict = self.data.remove_station("pb23")

        self.assertNotIn("pb23", new_mt_dict.keys())
        self.assertNotIn("pb23", new_data["station"])


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
