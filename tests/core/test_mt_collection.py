# -*- coding: utf-8 -*-
"""
Collection of MT stations

Created on Mon Jan 11 15:36:38 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import unittest

import numpy as np

from mtpy import MT, MTCollection

import mt_metadata


# =============================================================================
#
# =============================================================================
class TestMTCollection(unittest.TestCase):
    def setUp(self):
        fn_list = [
            value for key, value in mt_metadata.__dict__.items() if key.startswith("TF")
        ]

        self.mc = MTCollection()
        self.mc.initialize_collection("test_collection")

        self.mc.add_tf(fn_list)

    def test_filename(self):
        self.assertEqual(
            self.mc.working_directory.joinpath("test_collection.h5"),
            self.mc.mth5_filename,
        )

    def test_dataframe(self):
        with self.subTest("shape"):
            self.assertEqual(self.mc.dataframe.shape, (11, 14))
        with self.subTest("has columns"):
            self.assertListEqual(
                self.mc.dataframe.columns.to_list(),
                [
                    "station",
                    "survey",
                    "latitude",
                    "longitude",
                    "elevation",
                    "tf_id",
                    "units",
                    "has_impedance",
                    "has_tipper",
                    "has_covariance",
                    "period_min",
                    "period_max",
                    "hdf5_reference",
                    "station_hdf5_reference",
                ],
            )
        with self.subTest("stations"):
            self.assertListEqual(
                self.mc.dataframe.station.to_list(),
                [
                    "14-IEB0537A",
                    "NMX20",
                    "24",
                    "300",
                    "BP05",
                    "GEO858",
                    "Geoscience Australia",
                    "SAGE_2005",
                    "TEST01",
                    "YSW212abcdefghijkl",
                    "s08",
                ],
            )
        with self.subTest("latitude"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.latitude.to_numpy(),
                    np.array(
                        [
                            -22.82372222,
                            34.470528,
                            32.83331167,
                            34.727,
                            0.0,
                            22.69137833,
                            -23.05113333,
                            35.55,
                            -30.930285,
                            44.631,
                            -34.646,
                        ]
                    ),
                ).all(),
            )
        with self.subTest("longitude"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.longitude.to_numpy(),
                    np.array(
                        [
                            139.29469444,
                            -108.712288,
                            -107.08305667,
                            -115.735,
                            0.0,
                            139.70504,
                            139.46753333,
                            -106.28333333,
                            127.22923,
                            -110.44,
                            137.006,
                        ]
                    ),
                ).all(),
            )
        with self.subTest("tf_id"):
            self.assertListEqual(
                self.mc.dataframe.tf_id.to_list(),
                [
                    "14-IEB0537A",
                    "NMX20",
                    "24",
                    "300",
                    "BP05",
                    "GEO858",
                    "Geoscience Australia",
                    "SAGE_2005",
                    "TEST01",
                    "ysw212abcdefghijkl",
                    "s08",
                ],
            )
        with self.subTest("has_impedance"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.has_impedance.to_numpy(),
                    np.array(
                        [
                            True,
                            True,
                            True,
                            True,
                            True,
                            True,
                            True,
                            True,
                            True,
                            False,
                            True,
                        ],
                        dtype=bool,
                    ),
                ).all()
            )
        with self.subTest("has_tipper"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.has_tipper.to_numpy(),
                    np.array(
                        [
                            True,
                            True,
                            False,
                            True,
                            False,
                            True,
                            True,
                            True,
                            True,
                            True,
                            False,
                        ],
                        dtype=bool,
                    ),
                ).all()
            )
        with self.subTest("has_covariance"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.has_covariance.to_numpy(),
                    np.array(
                        [
                            True,
                            True,
                            False,
                            True,
                            False,
                            False,
                            True,
                            True,
                            False,
                            True,
                            False,
                        ],
                        dtype=bool,
                    ),
                ).all()
            )
        with self.subTest("period_min"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.period_min.to_numpy(),
                    np.array(
                        [
                            3.12500000e-03,
                            4.65455000e00,
                            9.76562500e-04,
                            1.16364000e00,
                            1.33333300e00,
                            5.15463918e-03,
                            1.00612732e-04,
                            4.19639110e-03,
                            1.21152720e-03,
                            1.81800000e-02,
                            7.93999902e-03,
                        ]
                    ),
                ).all()
            )
        with self.subTest("period_max"):
            self.assertTrue(
                np.isclose(
                    self.mc.dataframe.period_max.to_numpy(),
                    np.array(
                        [
                            2.94117647e03,
                            2.91271100e04,
                            4.26657565e01,
                            1.09226670e04,
                            6.45500000e01,
                            1.44927536e03,
                            1.02400262e00,
                            2.09731544e02,
                            1.21152749e03,
                            4.09600000e03,
                            2.73083324e03,
                        ]
                    ),
                ).all()
            )

    def tearDown(self):
        self.mc.mth5_collection.close_mth5()
        self.mc.mth5_filename.unlink()
