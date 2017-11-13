from unittest import TestCase
import numpy as np
import pytest

from mtpy.utils.gis_tools import ll_to_utm, utm_to_ll, project_point_ll2utm, project_point_utm2ll, transform_ll_to_utm, \
    transform_utm_to_ll


class TestGisTools(TestCase):
    def setUp(self):
        self.nref = 23  # 23, "WGS-84"
        self.lat = -34.3
        self.lon = 149.2

        self.zone = '55H'
        self.easting = 702466.54318647704
        self.northing = 9997808.8686107509

    def test_ll_to_utm(self):
        zone, easting, northing = ll_to_utm(self.nref, self.lat, self.lon)

        self.assertTrue(zone, self.zone)
        self.assertTrue(easting, self.easting)
        self.assertTrue(northing, self.northing)

    @pytest.mark.xfail
    def test_utm_to_ll(self):
        # test UTM to LL
        new_lat, new_lon = utm_to_ll(self.nref, self.northing, self.easting, self.zone)

        self.assertTrue(np.isclose(self.lat, new_lat))
        self.assertTrue(np.isclose(self.lon, new_lon))

    def test_project_point_ll2utm(self):
        zone, easting, northing = project_point_ll2utm(self.lat, self.lon)

        self.assertTrue(zone, self.zone)
        self.assertTrue(easting, self.easting)
        self.assertTrue(northing, self.northing)

    @pytest.mark.xfail
    def test_project_point_utm2ll(self):
        new_lat, new_lon = project_point_utm2ll(self.easting, self.northing, self.zone)

        self.assertTrue(np.isclose(self.lat, new_lat))
        self.assertTrue(np.isclose(self.lon, new_lon))

    def test_transform_ll_to_utm(self):
        zone, easting, northing = transform_ll_to_utm(self.lat, self.lon)

        self.assertTrue(zone, self.zone)
        self.assertTrue(easting, self.easting)
        self.assertTrue(northing, self.northing)

    @pytest.mark.xfail
    def test_transform_utm_to_ll(self):
        new_lat, new_lon = transform_utm_to_ll(self.easting, self.northing, self.zone)

        self.assertTrue(np.isclose(self.lat, new_lat))
        self.assertTrue(np.isclose(self.lon, new_lon))
