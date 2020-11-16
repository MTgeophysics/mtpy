"""
Test EDI module
"""

# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
import unittest


from mtpy.core.io import edi
from tests import TEST_MTPY_ROOT, make_temp_dir, EDI_DATA_DIR_BB

edi_path = Path(EDI_DATA_DIR_BB)
edi_list = list(edi_path.glob("*.edi"))
# =============================================================================
# Test EDI
# =============================================================================
class TestEDI(unittest.TestCase):
    """
    Testing mtpy.core.io.edi.Edi
    """

    def setUp(self):
        self.edi_obj = edi.Edi()
        
class TestHeader(unittest.TestCase):
    """
    Testing mtpy.core.io.edi.Header
    """
        
    def setUp(self):
        self.header = edi.Header()
        
    def test_fn(self):
        self.header.fn = edi_list[0]
        
        self.assertIsInstance(self.header.fn, Path)
        
    def test_location(self):
        
        lat_test_float = 40.0
        lon_test_float = -120.0
        lat_test_str = "40:00:00"
        lon_test_str = "-120:00:00"
        elev_test_float = 100.0
        elev_test_str = "100"
        
        self.header.lat = lat_test_float
        self.assertEqual(self.header.lat, lat_test_float)
        self.header.lat = lat_test_str
        self.assertEqual(self.header.lat, lat_test_float)
        
        self.header.lon = lon_test_float
        self.assertEqual(self.header.lon, lon_test_float)
        self.header.lon = lon_test_str
        self.assertEqual(self.header.lon, lon_test_float)
        
        self.header.elev= elev_test_float
        self.assertEqual(self.header.elev, elev_test_float)
        self.header.elev= elev_test_str
        self.assertEqual(self.header.elev, elev_test_float)
        
    def test_read_phx(self):
        self.header.fn = edi_list[0]
        self.header.read_header()
        
        lat_float = -30.939149166666667
        lon_float = 127.12636305555554
        
        self.assertAlmostEqual(self.header.lat, lat_float, 3)
        self.assertAlmostEqual(self.header.lon, lon_float, 3)

        
        

        
        
        
        


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
