# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 16:27:01 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import unittest

import pandas as pd
from mtpy.core import MTStations, MTLocation
from mtpy import MT

# =============================================================================


class TestMTStationGrid(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.east = 243900.352
        self.north = 4432069.056898517
        self.utm_epsg = 32611
        self.center = MTLocation(
            latitude=40.036594,
            longitude=-119.978167,
            utm_epsg=32611,
            model_east=245900.352,
            model_north=4436069.057,
        )
        dx = 1000
        dy = 2000
        count = 1
        mt_list = []
        for ii in range(5):
            for jj in range(5):
                mt_obj = MT(
                    east=(self.east + ii * dx),
                    north=(self.north + jj * dy),
                    utm_epsg=self.utm_epsg,
                    station=f"mt{count:02}",
                )
                count += 1
                mt_list.append(mt_obj)

        self.stations = MTStations(self.utm_epsg, mt_list=mt_list)

        self.station_locations = pd.DataFrame(
            {
                "survey": {
                    0: "0",
                    1: "0",
                    2: "0",
                    3: "0",
                    4: "0",
                    5: "0",
                    6: "0",
                    7: "0",
                    8: "0",
                    9: "0",
                    10: "0",
                    11: "0",
                    12: "0",
                    13: "0",
                    14: "0",
                    15: "0",
                    16: "0",
                    17: "0",
                    18: "0",
                    19: "0",
                    20: "0",
                    21: "0",
                    22: "0",
                    23: "0",
                    24: "0",
                },
                "station": {
                    0: "mt01",
                    1: "mt02",
                    2: "mt03",
                    3: "mt04",
                    4: "mt05",
                    5: "mt06",
                    6: "mt07",
                    7: "mt08",
                    8: "mt09",
                    9: "mt10",
                    10: "mt11",
                    11: "mt12",
                    12: "mt13",
                    13: "mt14",
                    14: "mt15",
                    15: "mt16",
                    16: "mt17",
                    17: "mt18",
                    18: "mt19",
                    19: "mt20",
                    20: "mt21",
                    21: "mt22",
                    22: "mt23",
                    23: "mt24",
                    24: "mt25",
                },
                "latitude": {
                    0: 39.99999999999098,
                    1: 40.01799481886136,
                    2: 40.035989568723195,
                    3: 40.05398424955559,
                    4: 40.0719788613376,
                    5: 40.00030250935993,
                    6: 40.0182975199738,
                    7: 40.036292461678464,
                    8: 40.05428733445315,
                    9: 40.07228213827704,
                    10: 40.00060383965571,
                    11: 40.01859904126739,
                    12: 40.03659417406883,
                    13: 40.05458923803936,
                    14: 40.072584233158295,
                    15: 40.00090399082604,
                    16: 40.01889938268981,
                    17: 40.036894705841924,
                    18: 40.05488996026181,
                    19: 40.0728851459289,
                    20: 40.001202962818844,
                    21: 40.019198544188924,
                    22: 40.03719405694554,
                    23: 40.055189501068234,
                    24: 40.073184876536565,
                },
                "longitude": {
                    0: -120.00000000034771,
                    1: -120.00078857255814,
                    2: -120.00157785431533,
                    3: -120.00236784633223,
                    4: -120.00315854932288,
                    5: -119.988300873993,
                    6: -119.98908638218944,
                    7: -119.9898725971845,
                    8: -119.9906595196884,
                    9: -119.99144715041241,
                    10: -119.97660157111596,
                    11: -119.9773840151211,
                    12: -119.97816716317647,
                    13: -119.97895101598958,
                    14: -119.97973557426889,
                    15: -119.96490209240297,
                    16: -119.9656814720402,
                    17: -119.96646155297907,
                    18: -119.96724233592424,
                    19: -119.96802382158145,
                    20: -119.95320243854052,
                    21: -119.95397875363396,
                    22: -119.9547557672801,
                    23: -119.9555334801809,
                    24: -119.95631189303934,
                },
                "elevation": {
                    0: 0.0,
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 0.0,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: 0.0,
                    21: 0.0,
                    22: 0.0,
                    23: 0.0,
                    24: 0.0,
                },
                "datum_epsg": {
                    0: "4326",
                    1: "4326",
                    2: "4326",
                    3: "4326",
                    4: "4326",
                    5: "4326",
                    6: "4326",
                    7: "4326",
                    8: "4326",
                    9: "4326",
                    10: "4326",
                    11: "4326",
                    12: "4326",
                    13: "4326",
                    14: "4326",
                    15: "4326",
                    16: "4326",
                    17: "4326",
                    18: "4326",
                    19: "4326",
                    20: "4326",
                    21: "4326",
                    22: "4326",
                    23: "4326",
                    24: "4326",
                },
                "east": {
                    0: 243900.352,
                    1: 243900.352,
                    2: 243900.352,
                    3: 243900.352,
                    4: 243900.352,
                    5: 244900.352,
                    6: 244900.352,
                    7: 244900.352,
                    8: 244900.352,
                    9: 244900.352,
                    10: 245900.352,
                    11: 245900.352,
                    12: 245900.352,
                    13: 245900.352,
                    14: 245900.352,
                    15: 246900.352,
                    16: 246900.352,
                    17: 246900.352,
                    18: 246900.352,
                    19: 246900.352,
                    20: 247900.352,
                    21: 247900.352,
                    22: 247900.352,
                    23: 247900.352,
                    24: 247900.352,
                },
                "north": {
                    0: 4432069.056898517,
                    1: 4434069.056898517,
                    2: 4436069.056898517,
                    3: 4438069.056898517,
                    4: 4440069.056898517,
                    5: 4432069.056898517,
                    6: 4434069.056898517,
                    7: 4436069.056898517,
                    8: 4438069.056898517,
                    9: 4440069.056898517,
                    10: 4432069.056898517,
                    11: 4434069.056898517,
                    12: 4436069.056898517,
                    13: 4438069.056898517,
                    14: 4440069.056898517,
                    15: 4432069.056898517,
                    16: 4434069.056898517,
                    17: 4436069.056898517,
                    18: 4438069.056898517,
                    19: 4440069.056898517,
                    20: 4432069.056898517,
                    21: 4434069.056898517,
                    22: 4436069.056898517,
                    23: 4438069.056898517,
                    24: 4440069.056898517,
                },
                "utm_epsg": {
                    0: "32611",
                    1: "32611",
                    2: "32611",
                    3: "32611",
                    4: "32611",
                    5: "32611",
                    6: "32611",
                    7: "32611",
                    8: "32611",
                    9: "32611",
                    10: "32611",
                    11: "32611",
                    12: "32611",
                    13: "32611",
                    14: "32611",
                    15: "32611",
                    16: "32611",
                    17: "32611",
                    18: "32611",
                    19: "32611",
                    20: "32611",
                    21: "32611",
                    22: "32611",
                    23: "32611",
                    24: "32611",
                },
                "model_east": {
                    0: -2000.0,
                    1: -2000.0,
                    2: -2000.0,
                    3: -2000.0,
                    4: -2000.0,
                    5: -1000.0,
                    6: -1000.0,
                    7: -1000.0,
                    8: -1000.0,
                    9: -1000.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 1000.0,
                    16: 1000.0,
                    17: 1000.0,
                    18: 1000.0,
                    19: 1000.0,
                    20: 2000.0,
                    21: 2000.0,
                    22: 2000.0,
                    23: 2000.0,
                    24: 2000.0,
                },
                "model_north": {
                    0: -4000.0,
                    1: -2000.0,
                    2: 0.0,
                    3: 2000.0,
                    4: 4000.0,
                    5: -4000.0,
                    6: -2000.0,
                    7: 0.0,
                    8: 2000.0,
                    9: 4000.0,
                    10: -4000.0,
                    11: -2000.0,
                    12: 0.0,
                    13: 2000.0,
                    14: 4000.0,
                    15: -4000.0,
                    16: -2000.0,
                    17: 0.0,
                    18: 2000.0,
                    19: 4000.0,
                    20: -4000.0,
                    21: -2000.0,
                    22: 0.0,
                    23: 2000.0,
                    24: 4000.0,
                },
                "model_elevation": {
                    0: 0.0,
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 0.0,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: 0.0,
                    21: 0.0,
                    22: 0.0,
                    23: 0.0,
                    24: 0.0,
                },
                "profile_offset": {
                    0: 0.0,
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 0.0,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: 0.0,
                    21: 0.0,
                    22: 0.0,
                    23: 0.0,
                    24: 0.0,
                },
            }
        )

        self.rotated_station_locations = pd.DataFrame(
            {
                "survey": {
                    0: "0",
                    1: "0",
                    2: "0",
                    3: "0",
                    4: "0",
                    5: "0",
                    6: "0",
                    7: "0",
                    8: "0",
                    9: "0",
                    10: "0",
                    11: "0",
                    12: "0",
                    13: "0",
                    14: "0",
                    15: "0",
                    16: "0",
                    17: "0",
                    18: "0",
                    19: "0",
                    20: "0",
                    21: "0",
                    22: "0",
                    23: "0",
                    24: "0",
                },
                "station": {
                    0: "mt01",
                    1: "mt02",
                    2: "mt03",
                    3: "mt04",
                    4: "mt05",
                    5: "mt06",
                    6: "mt07",
                    7: "mt08",
                    8: "mt09",
                    9: "mt10",
                    10: "mt11",
                    11: "mt12",
                    12: "mt13",
                    13: "mt14",
                    14: "mt15",
                    15: "mt16",
                    16: "mt17",
                    17: "mt18",
                    18: "mt19",
                    19: "mt20",
                    20: "mt21",
                    21: "mt22",
                    22: "mt23",
                    23: "mt24",
                    24: "mt25",
                },
                "latitude": {
                    0: 39.99999999999098,
                    1: 40.01799481886136,
                    2: 40.035989568723195,
                    3: 40.05398424955559,
                    4: 40.0719788613376,
                    5: 40.00030250935993,
                    6: 40.0182975199738,
                    7: 40.036292461678464,
                    8: 40.05428733445315,
                    9: 40.07228213827704,
                    10: 40.00060383965571,
                    11: 40.01859904126739,
                    12: 40.03659417406883,
                    13: 40.05458923803936,
                    14: 40.072584233158295,
                    15: 40.00090399082604,
                    16: 40.01889938268981,
                    17: 40.036894705841924,
                    18: 40.05488996026181,
                    19: 40.0728851459289,
                    20: 40.001202962818844,
                    21: 40.019198544188924,
                    22: 40.03719405694554,
                    23: 40.055189501068234,
                    24: 40.073184876536565,
                },
                "longitude": {
                    0: -120.00000000034771,
                    1: -120.00078857255814,
                    2: -120.00157785431533,
                    3: -120.00236784633223,
                    4: -120.00315854932288,
                    5: -119.988300873993,
                    6: -119.98908638218944,
                    7: -119.9898725971845,
                    8: -119.9906595196884,
                    9: -119.99144715041241,
                    10: -119.97660157111596,
                    11: -119.9773840151211,
                    12: -119.97816716317647,
                    13: -119.97895101598958,
                    14: -119.97973557426889,
                    15: -119.96490209240297,
                    16: -119.9656814720402,
                    17: -119.96646155297907,
                    18: -119.96724233592424,
                    19: -119.96802382158145,
                    20: -119.95320243854052,
                    21: -119.95397875363396,
                    22: -119.9547557672801,
                    23: -119.9555334801809,
                    24: -119.95631189303934,
                },
                "elevation": {
                    0: 0.0,
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 0.0,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: 0.0,
                    21: 0.0,
                    22: 0.0,
                    23: 0.0,
                    24: 0.0,
                },
                "datum_epsg": {
                    0: "4326",
                    1: "4326",
                    2: "4326",
                    3: "4326",
                    4: "4326",
                    5: "4326",
                    6: "4326",
                    7: "4326",
                    8: "4326",
                    9: "4326",
                    10: "4326",
                    11: "4326",
                    12: "4326",
                    13: "4326",
                    14: "4326",
                    15: "4326",
                    16: "4326",
                    17: "4326",
                    18: "4326",
                    19: "4326",
                    20: "4326",
                    21: "4326",
                    22: "4326",
                    23: "4326",
                    24: "4326",
                },
                "east": {
                    0: 243900.352,
                    1: 243900.352,
                    2: 243900.352,
                    3: 243900.352,
                    4: 243900.352,
                    5: 244900.352,
                    6: 244900.352,
                    7: 244900.352,
                    8: 244900.352,
                    9: 244900.352,
                    10: 245900.352,
                    11: 245900.352,
                    12: 245900.352,
                    13: 245900.352,
                    14: 245900.352,
                    15: 246900.352,
                    16: 246900.352,
                    17: 246900.352,
                    18: 246900.352,
                    19: 246900.352,
                    20: 247900.352,
                    21: 247900.352,
                    22: 247900.352,
                    23: 247900.352,
                    24: 247900.352,
                },
                "north": {
                    0: 4432069.056898517,
                    1: 4434069.056898517,
                    2: 4436069.056898517,
                    3: 4438069.056898517,
                    4: 4440069.056898517,
                    5: 4432069.056898517,
                    6: 4434069.056898517,
                    7: 4436069.056898517,
                    8: 4438069.056898517,
                    9: 4440069.056898517,
                    10: 4432069.056898517,
                    11: 4434069.056898517,
                    12: 4436069.056898517,
                    13: 4438069.056898517,
                    14: 4440069.056898517,
                    15: 4432069.056898517,
                    16: 4434069.056898517,
                    17: 4436069.056898517,
                    18: 4438069.056898517,
                    19: 4440069.056898517,
                    20: 4432069.056898517,
                    21: 4434069.056898517,
                    22: 4436069.056898517,
                    23: 4438069.056898517,
                    24: 4440069.056898517,
                },
                "utm_epsg": {
                    0: "32611",
                    1: "32611",
                    2: "32611",
                    3: "32611",
                    4: "32611",
                    5: "32611",
                    6: "32611",
                    7: "32611",
                    8: "32611",
                    9: "32611",
                    10: "32611",
                    11: "32611",
                    12: "32611",
                    13: "32611",
                    14: "32611",
                    15: "32611",
                    16: "32611",
                    17: "32611",
                    18: "32611",
                    19: "32611",
                    20: "32611",
                    21: "32611",
                    22: "32611",
                    23: "32611",
                    24: "32611",
                },
                "model_east": {
                    0: -4242.640687119285,
                    1: -2828.42712474619,
                    2: -1414.213562373095,
                    3: 0.0,
                    4: 1414.213562373095,
                    5: -3535.533905932738,
                    6: -2121.3203435596424,
                    7: -707.1067811865476,
                    8: 707.1067811865476,
                    9: 2121.3203435596424,
                    10: -2828.42712474619,
                    11: -1414.213562373095,
                    12: 0.0,
                    13: 1414.213562373095,
                    14: 2828.42712474619,
                    15: -2121.3203435596424,
                    16: -707.1067811865476,
                    17: 707.1067811865476,
                    18: 2121.3203435596424,
                    19: 3535.533905932738,
                    20: -1414.213562373095,
                    21: 0.0,
                    22: 1414.213562373095,
                    23: 2828.42712474619,
                    24: 4242.640687119285,
                },
                "model_north": {
                    0: -1414.213562373095,
                    1: 0.0,
                    2: 1414.213562373095,
                    3: 2828.42712474619,
                    4: 4242.640687119285,
                    5: -2121.3203435596424,
                    6: -707.1067811865476,
                    7: 707.1067811865476,
                    8: 2121.3203435596424,
                    9: 3535.533905932738,
                    10: -2828.42712474619,
                    11: -1414.213562373095,
                    12: 0.0,
                    13: 1414.213562373095,
                    14: 2828.42712474619,
                    15: -3535.533905932738,
                    16: -2121.3203435596424,
                    17: -707.1067811865476,
                    18: 707.1067811865476,
                    19: 2121.3203435596424,
                    20: -4242.640687119285,
                    21: -2828.42712474619,
                    22: -1414.213562373095,
                    23: 0.0,
                    24: 1414.213562373095,
                },
                "model_elevation": {
                    0: 0.0,
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 0.0,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: 0.0,
                    21: 0.0,
                    22: 0.0,
                    23: 0.0,
                    24: 0.0,
                },
                "profile_offset": {
                    0: 0.0,
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 0.0,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 0.0,
                    13: 0.0,
                    14: 0.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: 0.0,
                    21: 0.0,
                    22: 0.0,
                    23: 0.0,
                    24: 0.0,
                },
            }
        )

    def test_station_len(self):
        self.assertEqual(25, len(self.stations))

    def test_center_point(self):
        self.assertEqual(self.center, self.stations.center_point)

    def test_station_locations(self):
        self.assertTrue(
            (self.station_locations == self.stations.station_locations)
            .all()
            .all()
        )

    def test_copy(self):
        sc = self.stations.copy()
        self.assertEqual(self.stations, sc)

    def test_station_locations_rotated(self):
        s = self.stations.copy()
        s.rotate_stations(45)
        with self.subTest("rotated df"):
            self.assertTrue(
                (s.station_locations == self.rotated_station_locations)
                .all()
                .all()
            )
        with self.subTest("rotation_angle"):
            self.assertEqual(s.rotation_angle, 45)


class TestMTStationProfile(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.east = 243900.352
        self.north = 4432069.056898517
        self.utm_epsg = 32611
        self.center = MTLocation(
            latitude=42.212595,
            longitude=-120.078305,
            utm_epsg=32611,
            model_east=245900.352,
            model_north=4677969.409,
        )
        slope = 1
        count = 1
        dx = 1000
        mt_list = []
        for ii in range(5):
            x = self.east + ii * dx
            mt_obj = MT(
                east=x,
                north=slope * x + self.north,
                utm_epsg=self.utm_epsg,
                station=f"mt{count:02}",
            )
            count += 1
            mt_list.append(mt_obj)

        self.stations = MTStations(self.utm_epsg, mt_list=mt_list)

        self.profile_deg = (
            -120.10161927978938,
            42.19396005306167,
            -120.05497729522492,
            42.23123311383864,
            {"slope": 0.7991311074137044, "intercept": 138.17090007029887},
        )

        self.profile_m = (
            243900.352,
            4675969.408898517,
            247900.352,
            4679969.408898517,
            {"slope": 1.0, "intercept": 4432069.056898517},
        )

    def test_station_len(self):
        self.assertEqual(5, len(self.stations))

    def test_generate_profile_deg(self):
        self.assertTupleEqual(
            self.profile_deg, self.stations.generate_profile(units="deg")
        )

    def test_generate_profile_m(self):
        self.assertTupleEqual(
            self.profile_m, self.stations.generate_profile(units="m")
        )

    def test_center_point(self):
        self.assertEqual(self.center, self.stations.center_point)

    def test_extract_profile(self):
        d = self.stations._extract_profile(
            243900.352, 4675969.408898517, 247900.352, 4679969.408898517, 1000
        )
        d_names = [xx.station for xx in d]
        self.assertListEqual(["mt01", "mt02", "mt03", "mt04", "mt05"], d_names)


class TestMTStationMoreThanOneEPSG(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        m1 = MT(latitude=40, longitude=-120, utm_epsg=32611, station="mt01")
        m2 = MT(latitude=20, longitude=-100, utm_epsg=32613, station="mt02")
        m3 = MT(latitude=42, longitude=-118, utm_epsg=32611, station="mt03")

        self.stations = MTStations(32613, mt_list=[m1, m2, m3])

    def test_utm_epsg(self):
        self.assertEqual(self.stations.utm_epsg, 32611)


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
