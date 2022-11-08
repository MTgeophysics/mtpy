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

from mtpy.core.transfer_function.z import Z

# =============================================================================


class TestZInitialize(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.z = Z()

    def test_n_periods(self):
        self.assertEqual(0, self.z.n_periods)

    def test_is_empty(self):
        self.assertTrue(self.z._is_empty())

    def test_has_tf(self):
        self.assertFalse(self.z._has_tf())

    def test_has_tf_error(self):
        self.assertFalse(self.z._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertFalse(self.z._has_tf_model_error())

    def test_resistivity(self):
        self.assertEqual(None, self.z.resistivity)

    def test_phase(self):
        self.assertEqual(None, self.z.phase)


class TestZSetResPhase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.z = Z()
        self.resistivity = np.array([[[5.0, 100.0], [100.0, 5.0]]])
        self.resistivity_error = np.array([[1, 5], [5, 1]])
        self.phase_error = np.array([[0.5, 1], [1, 0.5]])
        self.phase = np.array([[[90.0, 45.0], [-135.0, -90.0]]])
        self.z.set_resistivity_phase(
            self.resistivity,
            self.phase,
            np.array([1]),
            res_err=self.resistivity_error,
            phase_err=self.phase_error,
        )

    def test_is_empty(self):
        self.assertFalse(self.z._is_empty())

    def test_has_tf(self):
        self.assertTrue(self.z._has_tf())

    def test_has_tf_error(self):
        self.assertTrue(self.z._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertFalse(self.z._has_tf_model_error())

    def test_resistivity(self):
        self.assertTrue(np.isclose(self.z.resistivity, self.resistivity).all())

    def test_phase(self):
        self.assertTrue(np.isclose(self.z.phase, self.phase).all())


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
