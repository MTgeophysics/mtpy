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
        self.utm_epsg = 32611
        self.utm_zone = "11N"
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

        with self.subTest("utm zone"):
            self.assertEqual(self.loc.utm_zone, self.utm_zone)


class TestMTLocationModelLocation(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.true_lat = 40.0
        self.true_lon = -120.0
        self.utm_epsg = 32611
        self.true_elevation = 1899.181396484

        self.loc = MTLocation(
            latitude=self.true_lat,
            longitude=self.true_lon,
            utm_epsg=self.utm_epsg,
        )
        self.center = MTLocation(
            latitude=self.true_lat,
            longitude=self.true_lon,
            utm_epsg=self.utm_epsg,
        )

        self.center.model_east = self.center.east
        self.center.model_north = self.center.north

        self.loc.compute_model_location(self.center)

    def test_model_location_east(self):
        self.assertEqual(self.loc.model_east, 0)

    def test_model_location_north(self):
        self.assertEqual(self.loc.model_north, 0)

    def test_get_elevation_from_national_map(self):
        self.loc.get_elevation_from_national_map()

        self.assertAlmostEqual(self.true_elevation, self.loc.elevation)

    def test_copy(self):
        loc_copy = self.loc.copy()
        self.assertEqual(loc_copy, self.loc)


class TestMTLocationModelLocation2(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.true_lat = 40.0
        self.true_lon = -120.0
        self.utm_epsg = 32611

        self.model_east = -431703.01876173366
        self.model_north = 224366.6894259695

        self.loc = MTLocation(
            latitude=self.true_lat,
            longitude=self.true_lon,
            utm_epsg=self.utm_epsg,
        )
        self.center = MTLocation(
            latitude=38,
            longitude=-115,
            utm_epsg=self.utm_epsg,
        )

        self.center.model_east = self.center.east
        self.center.model_north = self.center.north

        self.loc.compute_model_location(self.center)

    def test_model_location_east(self):
        self.assertAlmostEqual(self.loc.model_east, self.model_east)

    def test_model_location_north(self):
        self.assertAlmostEqual(self.loc.model_north, self.model_north)

    def test_project_onto_profile(self):
        self.loc.project_onto_profile_line(1, 240000)

        self.assertAlmostEqual(self.loc.profile_offset, 3136704.0501892385)


class TestMTLocationEqual(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.loc_01 = MTLocation(latitude=40, longitude=-120, utm_epsg=32611)
        self.loc_02 = MTLocation(
            east=243900.352029723, north=4432069.056898517, utm_epsg=32611
        )

    def test_equal(self):
        self.assertTrue(self.loc_01 == self.loc_02)

    def test_not_equal(self):
        loc_02 = MTLocation(
            east=244000.352029723, north=4432069.056898517, utm_epsg=32611
        )
        self.assertFalse(self.loc_01 == loc_02)

    def test_fail_wrong_object(self):
        def equals(value):
            return value == self.loc_01

        self.assertRaises(TypeError, equals, 10)


# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    unittest.main()
