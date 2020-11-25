from unittest import TestCase
import numpy as np
import pytest

from mtpy.utils import gis_tools


class TestGisTools(TestCase):
    def setUp(self):

        self.lat_hhmmss = "-34:17:57.99"
        self.lat_str = "-34.299442"
        self.lat_fail = "-34:29.9442"
        self.lat_d = -34.299442

        self.lon_hhmmss = "149:12:03.71"
        self.lon_str = "149.2010301"
        self.lon_fail = "149:12.0371"
        self.lon_d = 149.2010301

        self.elev_d = 1254.1
        self.elev_fail = "1200m"

        self.zone = "55H"
        self.zone_number = 55
        self.is_northern = False
        self.utm_letter = "H"
        self.zone_epsg = 32755
        self.easting = 702562.690286
        self.northing = 6202448.52785
        self.atol = 0.3  # tolerance of error
        self.from_epsg = 4326
        self.to_epsg = 28355

    def test_project_point_ll2utm(self):
        easting, northing, zone = gis_tools.project_point_ll2utm(self.lat_d, self.lon_d)

        if isinstance(zone, (np.bytes_, bytes)):
            zone = zone.decode("UTF-8")

        self.assertTrue(zone == self.zone)
        self.assertTrue(np.isclose(easting, self.easting))
        self.assertTrue(np.isclose(northing, self.northing))

    def test_project_point_utm2ll(self):
        new_lat, new_lon = gis_tools.project_point_utm2ll(
            self.easting, self.northing, self.zone
        )

        self.assertTrue(np.isclose(self.lat_d, new_lat))
        self.assertTrue(np.isclose(self.lon_d, new_lon))

        # testing with epsg
        new_lat, new_lon = gis_tools.project_point_utm2ll(
            self.easting, self.northing, utm_zone=self.zone, epsg=self.to_epsg
        )

        self.assertTrue(np.isclose(self.lat_d, new_lat))
        self.assertTrue(np.isclose(self.lon_d, new_lon))

    def test_convert_hhmmss(self):
        position_d = gis_tools.convert_position_str2float(self.lat_hhmmss)

        self.assertIsInstance(position_d, float)
        self.assertTrue(np.isclose(position_d, self.lat_d))

    def test_convert_str(self):
        position_d = gis_tools.convert_position_str2float(self.lat_str)

        self.assertIsInstance(position_d, float)
        self.assertTrue(np.isclose(position_d, self.lat_d))

    def test_convert_fail(self):
        with pytest.raises(gis_tools.GISError) as error:
            gis_tools.convert_position_str2float(self.position_fail)

    def test_convert_hhmmss(self):
        position_str = gis_tools.convert_position_float2str(self.lat_d)

        self.assertIsInstance(position_str, str)
        self.assertEqual(position_str, self.lat_hhmmss)

    def test_convert_fail(self):
        with pytest.raises(gis_tools.GISError) as error:
            gis_tools.convert_position_float2str(self.lat_fail)

    def test_assert_lat(self):
        lat_value = gis_tools.assert_lat_value(self.lat_str)

        self.assertIsInstance(lat_value, float)
        self.assertTrue(np.isclose(lat_value, self.lat_d))

    def test_assert_lon(self):
        lon_value = gis_tools.assert_lon_value(self.lon_str)

        self.assertIsInstance(lon_value, float)
        self.assertTrue(np.isclose(lon_value, self.lon_d))

    def test_get_utm_zone(self):
        zn, is_n, zone = gis_tools.get_utm_zone(self.lat_d, self.lon_d)

        self.assertEqual(zn, self.zone_number, "zone number")
        self.assertEqual(is_n, self.is_northern, "is northing")
        self.assertEqual(zone, self.zone, "utm zone")

    def test_utm_to_epsg(self):
        epsg_number = gis_tools.utm_zone_to_epsg(self.zone_number, self.is_northern)
        self.assertEqual(epsg_number, self.zone_epsg, "epsg number")

    def test_utm_letter_designation(self):
        utm_letter = gis_tools.utm_letter_designator(self.lat_hhmmss)

        self.assertEqual(utm_letter, self.utm_letter, "UTM letter")

    def test_validate_input_values(self):
        values = gis_tools.validate_input_values(self.lat_hhmmss, location_type="lat")

        self.assertIsInstance(values, np.ndarray)
        self.assertEqual(values.dtype.type, np.float64)

        values = gis_tools.validate_input_values(self.lon_hhmmss, location_type="lon")

        self.assertIsInstance(values, np.ndarray)
        self.assertEqual(values.dtype.type, np.float64)
