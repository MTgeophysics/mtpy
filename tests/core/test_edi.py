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

temp_dir = Path(make_temp_dir("EDI_TESTS"))
# =============================================================================
# Test EDI
# =============================================================================
class TestEDI(unittest.TestCase):
    """
    Testing mtpy.core.io.edi.Edi
    """

    def setUp(self):
        self.edi_obj = edi.Edi()
        
    def test_phoenix(self):
        self.edi_obj.fn = phoenix_fn
        
        self.assertEqual(self.edi_obj.station, "PHXTest01")
        self.assertAlmostEqual(self.edi_obj.lat, 10.123, 1)
        self.assertAlmostEqual(self.edi_obj.lon, 10.113, 1)
        self.assertAlmostEqual(self.edi_obj.elev, 2918, 1)
        self.assertEqual(self.edi_obj.Z.z.shape, (80, 2, 2))
        self.assertEqual(self.edi_obj.Tipper.tipper.shape, (80, 1, 2))
        self.assertEqual(self.edi_obj.Data.data_type_in, "spectra")
        
    def test_write_phoenix(self):
        mt_obj = edi.read_edi(phoenix_fn)
        edi_obj = edi.write_edi(mt_obj, fn=temp_dir.joinpath("phoenix_test.edi"))
        
        mt_obj = edi.read_edi(edi_obj.fn)
        
        self.assertEqual(edi_obj.station, "PHXTest01")
        self.assertEqual(mt_obj.ex_metadata, edi_obj.ex_metadata)
        self.assertEqual(mt_obj.ey_metadata, edi_obj.ey_metadata)
        self.assertEqual(mt_obj.hx_metadata, edi_obj.hx_metadata)
        self.assertEqual(mt_obj.hy_metadata, edi_obj.hy_metadata)
        self.assertEqual(mt_obj.hz_metadata, edi_obj.hz_metadata)
        
        
    def test_metronix(self):
        self.edi_obj.fn = metronix_fn
        
        self.assertEqual(self.edi_obj.station, "GEO858")
        self.assertAlmostEqual(self.edi_obj.lat, 22.691, 1)
        self.assertAlmostEqual(self.edi_obj.lon, 139.705, 1)
        self.assertAlmostEqual(self.edi_obj.elev, 181, 1)
        self.assertEqual(self.edi_obj.Z.z.shape, (73, 2, 2))
        self.assertEqual(self.edi_obj.Tipper.tipper.shape, (73, 1, 2))
        self.assertEqual(self.edi_obj.Data.data_type_in, "z")
        
    def test_write_metronix(self):
        mt_obj = edi.read_edi(metronix_fn)
        edi_obj = edi.write_edi(mt_obj, fn=temp_dir.joinpath("metronix_test.edi"))
        
        mt_obj = edi.read_edi(edi_obj.fn)
        
        self.assertEqual(edi_obj.station, "GEO858")
        self.assertEqual(mt_obj.ex_metadata, edi_obj.ex_metadata)
        self.assertEqual(mt_obj.ey_metadata, edi_obj.ey_metadata)
        self.assertEqual(mt_obj.hx_metadata, edi_obj.hx_metadata)
        self.assertEqual(mt_obj.hy_metadata, edi_obj.hy_metadata)
        self.assertEqual(mt_obj.hz_metadata, edi_obj.hz_metadata)
        
    def test_quantec(self):
        self.edi_obj.fn = quantec_fn
        
        self.assertEqual(self.edi_obj.station, "Geoscience Australia")
        self.assertAlmostEqual(self.edi_obj.lat, -23.051, 1)
        self.assertAlmostEqual(self.edi_obj.lon, 139.468, 1)
        self.assertAlmostEqual(self.edi_obj.elev, 122, 1)
        self.assertEqual(self.edi_obj.Z.z.shape, (41, 2, 2))
        self.assertEqual(self.edi_obj.Tipper.tipper.shape, (41, 1, 2))
        self.assertEqual(self.edi_obj.Data.data_type_in, "spectra")
        
    def test_write_quantec(self):
        mt_obj = edi.read_edi(quantec_fn)
        edi_obj = edi.write_edi(mt_obj, fn=temp_dir.joinpath("quantec_test.edi"))
        
        mt_obj = edi.read_edi(edi_obj.fn)
        
        self.assertEqual(edi_obj.station, "Geoscience Australia")
        self.assertEqual(mt_obj.ex_metadata, edi_obj.ex_metadata)
        self.assertEqual(mt_obj.ey_metadata, edi_obj.ey_metadata)
        self.assertEqual(mt_obj.hx_metadata, edi_obj.hx_metadata)
        self.assertEqual(mt_obj.hy_metadata, edi_obj.hy_metadata)
        self.assertEqual(mt_obj.hz_metadata, edi_obj.hz_metadata)


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
        
class TestDefine(unittest.TestCase):
    """ 
    Testing mtpy.core.io.edi.DefineMeasurement
    """
    
    def setUp(self):
        self.define = edi.DefineMeasurement()
        
    def test_phoenix(self):
        self.define.fn = phoenix_fn
        
        self.assertTrue(hasattr(self.define, "meas_ex"))
        self.assertTrue(hasattr(self.define, "meas_ey"))
        self.assertTrue(hasattr(self.define, "meas_hx"))
        self.assertTrue(hasattr(self.define, "meas_hy"))
        self.assertTrue(hasattr(self.define, "meas_hz"))
        self.assertTrue(hasattr(self.define, "meas_rrhx"))
        self.assertTrue(hasattr(self.define, "meas_rrhy"))
        
        self.assertEqual(self.define.meas_ex.chtype, "EX")
        self.assertEqual(self.define.meas_ey.chtype, "EY")
        self.assertEqual(self.define.meas_hx.chtype, "HX")
        self.assertEqual(self.define.meas_hy.chtype, "HY")
        self.assertEqual(self.define.meas_hz.chtype, "HZ")
        self.assertEqual(self.define.meas_rrhx.chtype, "RRHX")
        self.assertEqual(self.define.meas_rrhy.chtype, "RRHY")
        
        self.assertEqual(self.define.meas_ex.dipole_length, 50)
        self.assertEqual(self.define.meas_ex.azimuth, 0)
        
        self.assertAlmostEqual(self.define.meas_ey.dipole_length, 49.9, 1)
        self.assertAlmostEqual(self.define.meas_ey.azimuth, 116.67, 2)
        
        self.assertDictEqual(self.define.channel_ids,
                             {'EX': '114.011',
                             'EY': '115.011',
                             'HX': '111.011',
                             'HY': '112.011',
                             'HZ': '113.011',
                             'RRHX': '116.011',
                             'RRHY': '117.011'})
        
    def test_metronix(self):
        self.define.fn = metronix_fn
        
        self.assertTrue(hasattr(self.define, "meas_ex"))
        self.assertTrue(hasattr(self.define, "meas_ey"))
        self.assertTrue(hasattr(self.define, "meas_hx"))
        self.assertTrue(hasattr(self.define, "meas_hy"))
        self.assertTrue(hasattr(self.define, "meas_hz"))
        
        self.assertEqual(self.define.meas_ex.chtype, "EX")
        self.assertEqual(self.define.meas_ey.chtype, "EY")
        self.assertEqual(self.define.meas_hx.chtype, "HX")
        self.assertEqual(self.define.meas_hy.chtype, "HY")
        self.assertEqual(self.define.meas_hz.chtype, "HZ")
        
        
        self.assertEqual(self.define.meas_ex.dipole_length, 100)
        self.assertEqual(self.define.meas_ex.azimuth, 0)
        
        self.assertEqual(self.define.meas_ey.dipole_length, 100)
        self.assertEqual(self.define.meas_ey.azimuth, 90)
        
        self.assertDictEqual(self.define.channel_ids, 
                             {'EX': '1000.0001',
                             'EY': '1001.0001',
                             'HX': '1002.0001',
                             'HY': '1003.0001',
                             'HZ': '1004.0001'})
        
    def test_quantec(self):
        self.define.fn = quantec_fn
        
        self.assertTrue(hasattr(self.define, "meas_ex"))
        self.assertTrue(hasattr(self.define, "meas_ey"))
        self.assertTrue(hasattr(self.define, "meas_hx"))
        self.assertTrue(hasattr(self.define, "meas_hy"))
        self.assertTrue(hasattr(self.define, "meas_hz"))
        self.assertTrue(hasattr(self.define, "meas_rrhx"))
        self.assertTrue(hasattr(self.define, "meas_rrhy"))
        
        self.assertEqual(self.define.meas_ex.chtype, "EX")
        self.assertEqual(self.define.meas_ey.chtype, "EY")
        self.assertEqual(self.define.meas_hx.chtype, "HX")
        self.assertEqual(self.define.meas_hy.chtype, "HY")
        self.assertEqual(self.define.meas_hz.chtype, "HZ")
        self.assertEqual(self.define.meas_rrhx.chtype, "RRHX")
        self.assertEqual(self.define.meas_rrhy.chtype, "RRHY")
        
        self.assertEqual(self.define.meas_ex.dipole_length, 100)
        self.assertEqual(self.define.meas_ex.azimuth, 0)
        
        self.assertEqual(self.define.meas_ey.dipole_length, 100)
        self.assertEqual(self.define.meas_ey.azimuth, 90)
        
        self.assertDictEqual(self.define.channel_ids,
                             {'EX': '14.001',
                             'EY': '15.001',
                             'HX': '11.001',
                             'HY': '12.001',
                             'HZ': '13.001',
                             'RRHX': '11.001',
                             'RRHY': '12.001'})
        
class TestDataSection(unittest.TestCase):
    """
    Testing mtpy.core.io.edi.DataSection
    """
    def setUp(self):
        self.ds = edi.DataSection()
        
    def test_phoenix(self):
        self.ds.fn = phoenix_fn
        
        self.assertEqual(self.ds.data_type_in, "spectra")
        self.assertEqual(self.ds.nfreq, 80)
        self.assertEqual(self.ds.sectid, "PHXTest01")
        
    def test_metronix(self):
        self.ds.fn = metronix_fn
        
        self.assertEqual(self.ds.data_type_in, "z")
        self.assertEqual(self.ds.nfreq, 73)
        self.assertEqual(self.ds.sectid, "GEO858")
        
    def test_quantec(self):
        self.ds.fn = quantec_fn
        
        self.assertEqual(self.ds.data_type_in, "spectra")
        self.assertEqual(self.ds.nfreq, 41)
        self.assertEqual(self.ds.sectid, "IEA00184")
        

# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
