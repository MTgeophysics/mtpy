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
        self.mt_obj = mt.MT(self.j_fn)
        
    def test_location(self):
        self.assertEqual(self.j_obj.station_metadata.location.latitude,
                         self.mt_obj.latitude)
        self.assertEqual(self.j_obj.station_metadata.location.longitude,
                         self.mt_obj.longitude)
        self.assertEqual(self.j_obj.station_metadata.location.elevation,
                         self.mt_obj.elevation)
        
    def test_station(self):
        self.assertEqual(self.j_obj.station, self.mt_obj.station)
        

# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
    
        