# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:59:50 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

import numpy as np
from mtpy import MT
from mtpy.core.mt_dataframe import MTDataFrame

from mt_metadata import TF_EDI_CGG

# =============================================================================


class TestMT(unittest.TestCase):
    def setUp(self):
        self.mt = MT()

    def test_clone_empty(self):
        self.mt.station = "test_01"
        self.mt.survey = "big"
        self.mt.latitude = 10
        self.mt.longitude = 20
        new_mt = self.mt.clone_empty()

        for attr in ["survey", "station", "latitude", "longitude"]:
            with self.subTest(attr):
                self.assertEqual(getattr(new_mt, attr), getattr(self.mt, attr))

        with self.subTest("tf is empty"):
            self.assertFalse(new_mt.has_transfer_function())


class TestMTSetImpedance(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.z = np.array(
            [[0.1 - 0.1j, 10 + 10j], [-10 - 10j, -0.1 + 0.1j]]
        ).reshape((1, 2, 2))
        self.z_err = np.array([[0.1, 0.05], [0.05, 0.1]]).reshape((1, 2, 2))

        self.res = np.array([[[4.0e-03, 4.0e01], [4.0e01, 4.0e-03]]])
        self.res_err = np.array([[[0.002, 0.0005], [0.0005, 0.002]]])
        self.phase = np.array([[[-45.0, 45.0], [-135.0, 135.0]]])
        self.phase_err = np.array(
            [
                [
                    [2.65650512e01, 7.16197244e-04],
                    [7.16197244e-04, 2.65650512e01],
                ]
            ]
        )

        self.pt = np.array(
            [[[1.00020002, -0.020002], [-0.020002, 1.00020002]]]
        )

        self.mt = MT()
        self.mt.impedance = self.z
        self.mt.impedance_error = self.z_err
        self.mt.impedance_model_error = self.z_err

    def test_impedance(self):
        self.assertTrue((self.mt.impedance == self.z).all())

    def test_impedance_error(self):
        self.assertTrue(np.allclose(self.mt.impedance_error, self.z_err))

    def test_impedance_model_error(self):
        self.assertTrue(np.allclose(self.mt.impedance_model_error, self.z_err))

    def test_resistivity(self):
        self.assertTrue(np.allclose(self.mt.Z.resistivity, self.res))

    def test_resistivity_error(self):
        self.assertTrue(np.allclose(self.mt.Z.resistivity_error, self.res_err))

    def test_resistivity_model_error(self):
        self.assertTrue(
            np.allclose(self.mt.Z.resistivity_model_error, self.res_err)
        )

    def test_phase(self):
        self.assertTrue(np.allclose(self.mt.Z.phase, self.phase))

    def test_phase_error(self):
        self.assertTrue(np.allclose(self.mt.Z.phase_error, self.phase_err))

    def test_phase_model_error(self):
        self.assertTrue(
            np.allclose(self.mt.Z.phase_model_error, self.phase_err)
        )


class TestMT2DataFrame(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.m1 = MT(TF_EDI_CGG)
        self.m1.read_tf_file()

        self.sdf = self.m1.to_dataframe()
        self.mt_df = MTDataFrame()

    def test_station(self):
        self.assertEqual(self.sdf.station.unique()[0], "TEST01")

    def test_period(self):
        self.assertEqual(self.sdf.period.size, 73)

    def test_latitude(self):
        self.assertEqual(self.sdf.latitude.unique()[0], -30.930285)

    def test_longitude(self):
        self.assertEqual(self.sdf.longitude.unique()[0], 127.22923)

    def test_elevation(self):
        self.assertEqual(self.sdf.elevation.unique()[0], 175.27)

    def test_to_z(self):
        self.assertEqual(self.m1.Z, self.mt_df.to_z_object(self.sdf))

    def test_to_t(self):
        self.assertEqual(self.m1.Tipper, self.mt_df.to_t_object(self.sdf))


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
