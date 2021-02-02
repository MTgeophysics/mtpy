# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:54:38 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import unittest
from mtpy.core import mt
from mtpy.core.io.jfile import JFile
from tests import EDI_DATA_DIR_BB

# =============================================================================
# test j file
# =============================================================================
class TestJFile(unittest.TestCase):
    """
    test j file
    """

    def setUp(self):
        self.j_fn = list(EDI_DATA_DIR_BB.glob("*.j"))[0]
        self.j_obj = JFile(self.j_fn)
        
    def test_location(self):
        self.assertEqual(0.0, self.j_obj.latitude)
        self.assertEqual(0.0, self.j_obj.longitude)
        self.assertEqual(0.0, self.j_obj.elevation)
    
    def test_station(self):
        self.assertEqual("BP05", self.j_obj.station)
        
    def test_n_periods(self):
        # file has 2 null periods -999 so size is 12 not 14 like in the file. 
        self.assertEqual(12, self.j_obj.periods.size)
        
    def test_no_tipper(self):
        self.assertTrue(self.j_obj.Tipper.tipper.all() == 0.0)
        
        

class TestReadJFile(unittest.TestCase):
    """
    test reading j file into MT object
    """

    def setUp(self):
        self.j_fn = list(EDI_DATA_DIR_BB.glob("*.j"))[0]
        self.j_obj = JFile(self.j_fn)
        self.mt_obj = mt.MT(self.j_fn)

    def test_location(self):
        self.assertEqual(
            self.j_obj.station_metadata.location.latitude, self.mt_obj.latitude
        )
        self.assertEqual(
            self.j_obj.station_metadata.location.longitude, self.mt_obj.longitude
        )
        self.assertEqual(
            self.j_obj.station_metadata.location.elevation, self.mt_obj.elevation
        )

    def test_station(self):
        self.assertEqual(self.j_obj.station, self.mt_obj.station)
        
    def test_z(self):
        self.assertEqual(self.j_obj.Z, self.mt_obj.Z)
        
    def test_tipper(self):
        self.assertEqual(self.j_obj.Tipper, self.mt_obj.Tipper)

    def test_birrp_parameters(self):
        bp_list = [f"{k} = {v}" for k, v in self.j_obj.header_dict.items()]
        self.assertEqual(bp_list, 
                         self.mt_obj.station_metadata.transfer_function.processing_parameters)
        
        


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
