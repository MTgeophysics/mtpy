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
phoenix_fn = edi_list[-1]
quantec_fn = edi_list[2]
metronix_fn = edi_list[-2]
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

        self.header.elev = elev_test_float
        self.assertEqual(self.header.elev, elev_test_float)
        self.header.elev = elev_test_str
        self.assertEqual(self.header.elev, elev_test_float)

    def test_read_phx(self):
        self.header.fn = phoenix_fn
        self.header.read_header()

        lat_float = 10.122805555555557
        lon_float = 10.112722222222223

        self.assertAlmostEqual(self.header.lat, lat_float, 3)
        self.assertAlmostEqual(self.header.lon, lon_float, 3)

        self.assertEqual(self.header.acqdate, "2008-01-01")
        self.assertEqual(self.header.filedate, "2016-08-14")
        self.assertEqual(self.header.phoenix_edi, True)

    def test_read_metronix(self):
        self.header.fn = metronix_fn
        self.header.read_header()

        lat_float = 22.691378333333333
        lon_float = 139.70504

        self.assertAlmostEqual(self.header.lat, lat_float, 3)
        self.assertAlmostEqual(self.header.lon, lon_float, 3)

        self.assertEqual(self.header.acqdate, "2014-08-17")
        self.assertEqual(self.header.enddate, "2014-08-17")
        self.assertEqual(self.header.filedate, "2014-10-17")
        self.assertEqual(self.header.progdate, "2014-08-14")
        self.assertEqual(self.header.fileby, "Metronix")

    def test_read_quantec(self):
        self.header.fn = quantec_fn
        self.header.read_header()

        lat_float = -23.051133333333333
        lon_float = 139.46753333333334

        self.assertAlmostEqual(self.header.lat, lat_float, 3)
        self.assertAlmostEqual(self.header.lon, lon_float, 3)

        self.assertEqual(self.header.acqdate, "2014-11-15")
        self.assertEqual(self.header.enddate, "2014-11-15")
        self.assertEqual(self.header.filedate, "2014-11-17")
        self.assertEqual(self.header.progdate, "2012-10-10")
        self.assertEqual(self.header.fileby, "Quantec Geoscience")


class TestInformation(unittest.TestCase):
    """
    Testing mtpy.core.io.edi.Information
    """

    def setUp(self):
        self.info = edi.Information()

    def test_phoenix(self):
        self.info.fn = phoenix_fn

        self.assertNotEqual(self.info.info_dict, {})
        self.assertNotEqual(self.info.info_list, [])

    def test_metronix(self):
        self.info.fn = metronix_fn

        self.assertNotEqual(self.info.info_dict, {})
        self.assertNotEqual(self.info.info_list, [])

    def test_quantec(self):
        self.info.fn = quantec_fn

        self.assertEqual(self.info.info_dict, {})
        self.assertEqual(self.info.info_list, [])


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
