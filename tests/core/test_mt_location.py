# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:17:37 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

from mtpy.core.mt_location import MTLocation

# =============================================================================


class TestMTLocation(unittest.TestCase):
    def setUp(self):
        self.loc = MTLocation()
        self.true_lat = 40.0
        self.true_lon = -120.0
        self.true_elevation = 1200.0
        self.utm_epsg = 32611
        self.true_east = 243900.352029723
        self.true_north = 4432069.056898517

    def test_set_latitude_from_string_time(self):
        self.loc.latitude = "40:00:00"
        self.assertEqual(self.loc.latitude, self.true_lat)

    def test_set_latitude_from_string_decimal(self):
        self.loc.latitude = "40.0"
        self.assertEqual(self.loc.latitude, self.true_lat)

    def test_set_latitude_from_int(self):
        self.loc.latitude = int(self.true_lat)
        self.assertIsInstance(self.loc.latitude, float)

    def test_set_latitude_fail(self):

        with self.subTest("bad input"):
            self.assertRaises(ValueError, MTLocation, **{"latitude": "ten"})

        with self.subTest("Too large"):
            self.assertRaises(ValueError, MTLocation, **{"latitude": "100"})

    def test_set_longitude_from_string_time(self):
        self.loc.longitude = "-120:00:00"
        self.assertEqual(self.loc.longitude, self.true_lon)

    def test_set_longitude_from_string_decimal(self):
        self.loc.longitude = "-120.0"
        self.assertEqual(self.loc.longitude, self.true_lon)

    def test_set_longitude_from_int(self):
        self.loc.longitude = int(self.true_lon)
        self.assertIsInstance(self.loc.longitude, float)

    def test_set_longitude_fail(self):
        def set_lon(value):
            self.loc.longitude = value

        with self.subTest("bad input"):
            self.assertRaises(ValueError, MTLocation, **{"longitude": "ten"})

        with self.subTest("Too large"):
            self.assertRaises(ValueError, MTLocation, **{"longitude": "400"})

    def test_set_utm_zone(self):
        self.loc.latitude = self.true_lat
        self.loc.longitude = self.true_lon
        self.loc.utm_epsg = self.utm_epsg

        with self.subTest("east"):
            self.assertAlmostEqual(self.loc.east, self.true_east)
        with self.subTest("north"):
            self.assertAlmostEqual(self.loc.north, self.true_north)

    def test_set_utm_coordinates_fail(self):
        self.assertRaises(
            ValueError,
            MTLocation,
            **{"east": self.true_east, "north": self.true_north}
        )

    def test_set_utm_coordinates(self):
        self.loc = MTLocation(
            east=self.true_east, north=self.true_north, utm_epsg=self.utm_epsg
        )

        with self.subTest("latitude"):
            self.assertAlmostEqual(self.loc.latitude, self.true_lat)
        with self.subTest("longitude"):
            self.assertAlmostEqual(self.loc.longitude, self.true_lon)


# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    unittest.main()
