from unittest import TestCase
import numpy as np
import pytest

import gis_tools


class TestGisTools(TestCase):
    def setUp(self):
        self.lat = -34.299442
        self.lon = 149.201031

        self.zone = '55H'
        self.easting = 702562.773
        self.northing = 6202448.526
        self.atol = 0.3  # tolerance of error
        self.from_epsg = 4326
        self.to_epsg = 28355

    def test_project_point_ll2utm(self):
        easting, northing, zone = gis_tools.project_point_ll2utm(self.lat,
                                                                 self.lon)

        print((zone, easting, northing))


        if isinstance(zone, np.bytes_) or (zone, np.unicode_):
            zone = zone.decode('UTF-8')

        self.assertTrue(zone == self.zone)
        self.assertTrue(np.isclose(easting, self.easting))
        self.assertTrue(np.isclose(northing, self.northing))

    def test_project_point_utm2ll(self):
        new_lat, new_lon = gis_tools.project_point_utm2ll(self.easting,
                                                          self.northing,
                                                          self.zone)

        print((new_lat, new_lon))

        self.assertTrue(np.isclose(self.lat, new_lat))
        self.assertTrue(np.isclose(self.lon, new_lon))

        # testing with epsg
        new_lat, new_lon = gis_tools.project_point_utm2ll(self.easting, 
                                                          self.northing, 
                                                          utm_zone=self.zone, 
                                                          epsg=self.to_epsg)

        print((new_lat, new_lon))

        self.assertTrue(np.isclose(self.lat, new_lat))
        self.assertTrue(np.isclose(self.lon, new_lon))
