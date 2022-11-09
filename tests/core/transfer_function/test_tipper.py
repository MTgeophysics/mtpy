# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 13:04:38 2022

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import unittest
import numpy as np

from mtpy.core.transfer_function.tipper import Tipper

# =============================================================================


class TestTipperInitialize(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tipper = Tipper()

    def test_n_periods(self):
        self.assertEqual(0, self.tipper.n_periods)

    def test_is_empty(self):
        self.assertTrue(self.tipper._is_empty())

    def test_has_tf(self):
        self.assertFalse(self.tipper._has_tf())

    def test_has_tf_error(self):
        self.assertFalse(self.tipper._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertFalse(self.tipper._has_tf_model_error())

    def test_mag_real(self):
        self.assertEqual(None, self.tipper.mag_real)

    def test_mag_imag(self):
        self.assertEqual(None, self.tipper.mag_imag)

    def test_angle_real(self):
        self.assertEqual(None, self.tipper.angle_real)

    def test_angle_imag(self):
        self.assertEqual(None, self.tipper.angle_imag)

    def test_mag_error(self):
        self.assertEqual(None, self.tipper.mag_error)

    def test_mag_model_error(self):
        self.assertEqual(None, self.tipper.mag_model_error)

    def test_angle_error(self):
        self.assertEqual(None, self.tipper.angle_error)

    def test_angle_model_error(self):
        self.assertEqual(None, self.tipper.angle_model_error)

    def test_amplitude(self):
        self.assertEqual(None, self.tipper.amplitude)

    def test_phase(self):
        self.assertEqual(None, self.tipper.phase)

    def test_amplitude_error(self):
        self.assertEqual(None, self.tipper.amplitude_error)

    def test_phase_error(self):
        self.assertEqual(None, self.tipper.phase_error)

    def test_amplitude_model_error(self):
        self.assertEqual(None, self.tipper.amplitude_model_error)

    def test_phase_model_error(self):
        self.assertEqual(None, self.tipper.phase_model_error)


class TestSetTipper(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.t = np.ones((1, 1, 2)) + 0.25j * np.ones((1, 1, 2))
        self.t_error = np.ones((1, 1, 2)) * 0.01
        self.t_model_error = np.ones((1, 1, 2)) * 0.03
        self.tipper = Tipper(
            tipper=self.t,
            tipper_error=self.t_error,
            tipper_model_error=self.t_model_error,
            frequency=np.array([1]),
        )

    def test_is_empty(self):
        self.assertFalse(self.tipper._is_empty())

    def test_has_tf(self):
        self.assertTrue(self.tipper._has_tf())

    def test_has_tf_error(self):
        self.assertTrue(self.tipper._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertTrue(self.tipper._has_tf_model_error())

    def test_tipper(self):
        self.assertTrue(np.isclose(self.tipper.tipper, self.t).all())

    def test_tipper_error(self):
        self.assertTrue(
            np.isclose(self.tipper.tipper_error, self.t_error).all()
        )

    def test_tipper_model_error(self):
        self.assertTrue(
            np.isclose(
                self.tipper.tipper_model_error, self.t_model_error
            ).all()
        )

    def test_mag_real(self):
        self.assertTrue(
            np.isclose(self.tipper.mag_real, np.array([1.41421356])).all()
        )

    def test_mag_imag(self):
        self.assertTrue(
            np.isclose(self.tipper.mag_imag, np.array([0.35355339])).all()
        )

    def test_angle_real(self):
        self.assertTrue(
            np.isclose(self.tipper.angle_real, np.array([45.0])).all()
        )

    def test_angle_imag(self):
        self.assertTrue(
            np.isclose(self.tipper.angle_imag, np.array([45.0])).all()
        )

    def test_mag_error(self):
        self.assertTrue(
            np.isclose(self.tipper.mag_error, np.array([0.01414214])).all()
        )

    def test_angle_error(self):
        self.assertTrue(
            np.isclose(self.tipper.angle_error, np.array([0])).all()
        )

    def test_mag_model_error(self):
        self.assertTrue(
            np.isclose(
                self.tipper.mag_model_error, np.array([0.04242641])
            ).all()
        )

    def test_angle_model_error(self):
        self.assertTrue(
            np.isclose(self.tipper.angle_model_error, np.array([0])).all()
        )

    def test_amplitude(self):
        self.assertTrue(
            np.isclose(
                self.tipper.amplitude, np.array([[[1.03077641, 1.03077641]]])
            ).all()
        )

    def test_amplitude_error(self):
        self.assertTrue(
            np.isclose(
                self.tipper.amplitude_error,
                np.array([[[0.01212648, 0.01212648]]]),
            ).all()
        )

    def test_amplitude_model_error(self):
        self.assertTrue(
            np.isclose(
                self.tipper.amplitude_model_error,
                np.array([[[0.03637218, 0.03637218]]]),
            ).all()
        )

    def test_phase(self):
        self.assertTrue(
            np.isclose(
                self.tipper.phase, np.array([[[14.03624347, 14.03624347]]])
            ).all()
        )

    def test_phase_error(self):
        self.assertTrue(
            np.isclose(
                self.tipper.phase_error,
                np.array([[[0.67407048, 0.67407048]]]),
            ).all()
        )

    def test_phase_model_error(self):
        self.assertTrue(
            np.isclose(
                self.tipper.phase_model_error,
                np.array([[[2.02226993, 2.02226993]]]),
            ).all()
        )

    def test_rotate(self):
        b = self.tipper.rotate(45)

        with self.subTest("angle_real"):
            self.assertTrue(np.isclose(b.angle_real, np.array([90])).all())
        with self.subTest("angle_imag"):
            self.assertTrue(np.isclose(b.angle_imag, np.array([90])).all())
        with self.subTest("mag_error"):
            self.assertTrue(np.isclose(b.mag_error, np.array([0.02])).all())
        with self.subTest("mag_model_error"):
            self.assertTrue(
                np.isclose(b.mag_model_error, np.array([0.06])).all()
            )


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
