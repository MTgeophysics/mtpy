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
from mtpy.mtpy_globals import EDI_DATA_DIR, EDI_DATA_DIR2


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
        station_list = [fn.stem for fn in EDI_DATA_DIR.glob("*.edi")]

        self.assertListEqual(station_list, self.data.station_locations.station.tolist())

    def test_period_range(self):
        self.assertEqual(20, self.data.period_list.size)

    def test_add_station(self):
        new_data, new_mt_dict = self.data.add_station(
            fn=list(EDI_DATA_DIR2.glob("*.edi"))[0:2]
        )


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
