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

from mtpy.core.transfer_function.pt import PhaseTensor

# =============================================================================


class TestPTInitialize(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pt = PhaseTensor()

    def test_n_periods(self):
        self.assertEqual(0, self.pt.n_periods)

    def test_is_empty(self):
        self.assertTrue(self.pt._is_empty())

    def test_has_tf(self):
        self.assertFalse(self.pt._has_tf())

    def test_has_tf_error(self):
        self.assertFalse(self.pt._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertFalse(self.pt._has_tf_model_error())

    def test_empty_properties(self):
        for attr in [
            "pt",
            "phimin",
            "phimax",
            "alpha",
            "beta",
            "skew",
            "trace",
            "azimuth",
            "ellipticity",
            "eccentricity",
        ]:
            with self.subTest(f"{attr}"):
                self.assertEqual(None, getattr(self.pt, attr))
            with self.subTest(f"{attr}_error"):
                self.assertEqual(None, getattr(self.pt, f"{attr}_error"))
            with self.subTest(f"{attr}_model_error"):
                self.assertEqual(None, getattr(self.pt, f"{attr}_model_error"))


class TestZSetResPhase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        z = np.array([[0, 1 + 1j], [-1 - 1j, 0]])
        z_error = np.array([[0.1, 0.05], [0.05, 0.1]])
        self.pt = PhaseTensor(z=z, z_error=z_error, z_model_error=z_error)

    def test_is_empty(self):
        self.assertFalse(self.pt._is_empty())

    def test_has_tf(self):
        self.assertTrue(self.pt._has_tf())

    def test_has_tf_error(self):
        self.assertTrue(self.pt._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertTrue(self.pt._has_tf_model_error())

    def test_pt(self):
        self.assertTrue(
            np.isclose(self.pt.pt, np.array([[[1.0, 0.0], [0.0, 1.0]]])).all()
        )

    def test_pt_error(self):
        self.assertTrue(
            np.isclose(
                self.pt.pt_error, np.array([[[0.1, 0.2], [0.2, 0.1]]])
            ).all()
        )

    def test_pt_model_error(self):
        self.assertTrue(
            np.isclose(
                self.pt.pt_model_error, np.array([[[0.1, 0.2], [0.2, 0.1]]])
            ).all()
        )

    def test_det(self):
        self.assertTrue(np.isclose(self.pt.det, np.array([1.0])).all())

    def test_phimax(self):
        self.assertTrue(np.isclose(self.pt.phimax, np.array([45.0])).all())

    def test_phimin(self):
        self.assertTrue(np.isclose(self.pt.phimin, np.array([45.0])).all())

    def test_azimuth(self):
        self.assertTrue(np.isclose(self.pt.azimuth, np.array([0.0])).all())

    def test_beta(self):
        self.assertTrue(np.isclose(self.pt.beta, np.array([0.0])).all())

    def test_skew(self):
        self.assertTrue(np.isclose(self.pt.skew, np.array([0.0])).all())

    def test_ellipticity(self):
        self.assertTrue(np.isclose(self.pt.ellipticity, np.array([0.0])).all())

    def test_eccentricity(self):
        self.assertTrue(np.isclose(self.pt.eccentricity, np.array([0.0])).all())

    def test_det_error(self):
        self.assertTrue(np.isclose(self.pt.det_error, np.array([0.2])).all())

    def test_phimax_error(self):
        self.assertTrue(np.all(np.isnan(self.pt.phimax_error)))

    def test_phimin_error(self):
        self.assertTrue(np.all(np.isnan(self.pt.phimin_error)))

    def test_azimuth_error(self):
        self.assertTrue(np.all(np.isnan(self.pt.azimuth_error)))

    def test_beta_error(self):
        self.assertTrue(
            np.isclose(self.pt.beta_error, np.array([4.05142342])).all()
        )

    def test_skew_error(self):
        self.assertTrue(
            np.isclose(self.pt.skew_error, np.array([4.05142342])).all()
        )

    def test_ellipticity_error(self):
        self.assertTrue(np.all(np.isnan(self.pt.ellipticity_error)))

    def test_eccentricity_error(self):
        self.assertTrue(np.all(np.isnan(self.pt.eccentricity_error)))


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
