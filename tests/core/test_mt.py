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
    @classmethod
    def setUpClass(self):
        self.mt = MT()
        self.mt.station = "test_01"
        self.mt.survey = "big"
        self.mt.latitude = 10
        self.mt.longitude = 20

    def test_clone_empty(self):
        new_mt = self.mt.clone_empty()

        for attr in ["survey", "station", "latitude", "longitude"]:
            with self.subTest(attr):
                self.assertEqual(getattr(new_mt, attr), getattr(self.mt, attr))

        with self.subTest("tf is empty"):
            self.assertFalse(new_mt.has_transfer_function())

    def test_copy(self):
        mt_copy = self.mt.copy()

        self.assertEqual(self.mt, mt_copy)


class TestMTFromKWARGS(unittest.TestCase):
    def setUp(self):
        self.mt = MT(east=243900.352, north=4432069.056898517, utm_epsg=32611)

    def test_latitude(self):
        self.assertAlmostEqual(self.mt.latitude, 40)

    def test_longitude(self):
        self.assertAlmostEqual(self.mt.longitude, -120)


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
        self.pt_error = np.array(
            [[[0.01040308, 0.02020604], [0.02020604, 0.01040308]]]
        )
        self.pt_azimuth = np.array([315.0])
        self.pt_azimuth_error = np.array([3.30832308])
        self.pt_skew = np.array([0])
        self.pt_skew_error = np.array([0.40923428])

        self.mt = MT()
        self.mt.station = "mt001"
        self.mt.impedance = self.z
        self.mt.impedance_error = self.z_err
        self.mt.impedance_model_error = self.z_err

    def test_period(self):
        self.assertTrue((np.array([1]) == self.mt.period).all())

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

    def test_phase_tensor(self):
        self.assertTrue(np.allclose(self.pt, self.mt.pt.pt))

    def test_phase_tensor_error(self):
        self.assertTrue(np.allclose(self.pt_error, self.mt.pt.pt_error))

    def test_phase_tensor_model_error(self):
        self.assertTrue(np.allclose(self.pt_error, self.mt.pt.pt_model_error))

    def test_phase_tensor_azimuth(self):
        self.assertTrue(np.allclose(self.pt_azimuth, self.mt.pt.azimuth))

    def test_phase_tensor_azimuth_error(self):
        self.assertTrue(
            np.allclose(self.pt_azimuth_error, self.mt.pt.azimuth_error)
        )

    def test_phase_tensor_azimuth_model_error(self):
        self.assertTrue(
            np.allclose(self.pt_azimuth_error, self.mt.pt.azimuth_model_error)
        )

    def test_phase_tensor_skew(self):
        self.assertTrue(np.allclose(self.pt_skew, self.mt.pt.skew))

    def test_phase_tensor_skew_error(self):
        self.assertTrue(np.allclose(self.pt_skew_error, self.mt.pt.skew_error))

    def test_phase_tensor_skew_model_error(self):
        self.assertTrue(
            np.allclose(self.pt_skew_error, self.mt.pt.skew_model_error)
        )

    def test_remove_static_shift(self):
        new_mt = self.mt.remove_static_shift(ss_x=0.5, ss_y=1.5, inplace=False)

        self.assertTrue(
            np.allclose(
                (self.mt.impedance.data / new_mt.impedance.data) ** 2,
                np.array(
                    [[[0.5 + 0.0j, 0.5 + 0.0j], [1.5 - 0.0j, 1.5 - 0.0j]]]
                ),
            )
        )

    def test_remove_distortion(self):
        new_mt = self.mt.remove_distortion()

        self.assertTrue(
            np.all(
                np.isclose(
                    new_mt.Z.z,
                    np.array(
                        [
                            [
                                [
                                    0.099995 - 0.099995j,
                                    9.99949999 + 9.99949999j,
                                ],
                                [
                                    -9.99949999 - 9.99949999j,
                                    -0.099995 + 0.099995j,
                                ],
                            ]
                        ]
                    ),
                )
            )
        )

    def test_interpolate_fail_bad_f_type(self):
        self.assertRaises(
            ValueError, self.mt.interpolate, [0, 1], f_type="wrong"
        )

    def test_interpolate_fail_bad_periods(self):
        self.assertRaises(ValueError, self.mt.interpolate, [0.1, 2])


class TestMTComputeModelError(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.z = np.array(
            [[0.1 - 0.1j, 10 + 10j], [-10 - 10j, -0.1 + 0.1j]]
        ).reshape((1, 2, 2))
        self.z_err = np.array([[0.1, 0.05], [0.05, 0.1]]).reshape((1, 2, 2))

        self.mt = MT()
        self.mt.impedance = self.z
        self.mt.impedance_error = self.z_err

    def test_compute_model_error(self):
        err = np.array([[[0.70710678, 0.70710678], [0.70710678, 0.70710678]]])

        self.mt.compute_model_z_errors()

        self.assertTrue(np.allclose(self.mt.impedance_model_error.data, err))


class TestSetTipper(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.t = np.array([[[0.25 - 0.2j, 0.25 + 0.2j]]])
        self.t_err = np.array([[[0.02, 0.03]]])

        self.mt = MT()
        self.mt.tipper = self.t
        self.mt.tipper_error = self.t_err

    def test_tipper(self):
        self.assertTrue(np.allclose(self.mt.tipper.data, self.t))

    def test_tipper_error(self):
        self.assertTrue(np.allclose(self.mt.tipper_error.data, self.t_err))

    def test_tipper_model_error(self):
        err = np.array([[[0.02, 0.03]]])
        self.mt.compute_model_t_errors(
            error_type="absolute", error_value=0.02, floor=True
        )
        self.assertTrue(np.allclose(self.mt.tipper_model_error.data, err))


class TestMT2DataFrame(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.m1 = MT(TF_EDI_CGG)
        self.m1.read()

        self.mt_df = self.m1.to_dataframe()

    def test_station(self):
        self.assertEqual(self.mt_df.station, "TEST01")

    def test_period(self):
        self.assertEqual(self.mt_df.period.size, 73)

    def test_latitude(self):
        self.assertEqual(self.mt_df.latitude, -30.930285)

    def test_longitude(self):
        self.assertEqual(self.mt_df.longitude, 127.22923)

    def test_elevation(self):
        self.assertEqual(self.mt_df.elevation, 175.27)

    def test_to_z(self):
        self.assertEqual(self.m1.Z, self.mt_df.to_z_object())

    def test_to_t(self):
        self.assertEqual(self.m1.Tipper, self.mt_df.to_t_object())

    def test_isinstance_mt_dataframe(self):
        self.assertIsInstance(self.mt_df, MTDataFrame)

    def test_from_dataframe(self):
        m2 = MT()
        m2.from_dataframe(self.mt_df)

        for key in [
            "station",
            "latitude",
            "longitude",
            "elevation",
            "east",
            "north",
            "utm_epsg",
            "model_north",
            "model_east",
            "model_elevation",
        ]:
            with self.subTest(key):
                self.assertTrue(getattr(m2, key) == getattr(self.m1, key))

        with self.subTest("dataset"):
            self.assertEqual(self.m1._transfer_function, m2._transfer_function)

    def test_from_dataframe_fail(self):
        self.assertRaises(TypeError, self.m1.from_dataframe, "a")


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
