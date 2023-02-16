# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:09:59 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import unittest

import numpy as np
from mtpy.modeling.errors import ModelErrors

# =============================================================================


class TestErrorEstimation(unittest.TestCase):
    def setUp(self):

        self.data = np.array(
            [
                [[0.1 - 0.1j, 10 + 10j], [-0.3 + 0.3j, -20 - 20j]],
                [[0.2 - 0.2j, 5 + 5j], [-0.6 + 0.6j, -40 - 40j]],
            ]
        )
        self.data_error = np.array(
            [[[0.1, 0.02], [0.07, 0.25]], [[0.05, 0.1], [0.03, 0.025]]]
        )
        self.m = ModelErrors(self.data, self.data_error)

    def test_bad_shape(self):

        self.assertRaises(
            ValueError, self.m.validate_array_shape, np.random.rand(3, 3, 3)
        )

    def test_bad_shape_error(self):

        self.assertRaises(
            ValueError, self.m.validate_array_shape, np.random.rand(3, 3, 3)
        )

    def test_set_error_type_fail(self):
        def set_error_type(item):
            self.m.error_type = item

        self.assertRaises(NotImplementedError, set_error_type, "fail")

    def test_set_error_mode_fail(self):
        def set_mode(item):
            self.m.mode = item

        self.assertRaises(NotImplementedError, set_mode, "fail")

    def test_set_floor_fail(self):
        def set_floor(item):
            self.m.floor = item

        self.assertRaises(ValueError, set_floor, "fail")

    def test_set_error_value_percent(self):
        self.m.error_value = 5
        self.assertEqual(0.05, self.m.error_value)

    def test_compute_percent_error(self):
        err = self.m.compute_error(
            error_type="percent", error_value=5, floor=False
        )
        self.assertTrue(np.allclose(np.abs(self.data) * 0.05, err))

    def test_compute_percent_error_floor(self):
        err = self.m.compute_error(
            error_type="percent", error_value=5, floor=True
        )
        self.assertFalse(np.allclose(np.abs(self.data) * 0.05, err))

    def test_compute_off_diagonal_mean_error(self):
        err = self.m.compute_error(error_type="arithmetic_mean", floor=False)

        est = np.array(
            [
                [[0.35371245, 0.35371245], [0.35371245, 0.35371245]],
                [[0.17804494, 0.17804494], [0.17804494, 0.17804494]],
            ]
        )

        self.assertTrue(np.allclose(err, est))

    def test_compute_off_diagonal_mean_error_floor(self):
        err = self.m.compute_error(error_type="arithmetic_mean", floor=True)

        est = np.array(
            [
                [[0.35371245, 0.35371245], [0.35371245, 0.35371245]],
                [[0.17804494, 0.17804494], [0.17804494, 0.17804494]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_median(self):
        err = self.m.compute_error(error_type="median", floor=False)

        est = np.array(
            [
                [[0.00707107, 0.00707107], [0.00707107, 0.00707107]],
                [[0.01414214, 0.01414214], [0.01414214, 0.01414214]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_median_floor(self):
        err = self.m.compute_error(error_type="median", floor=True)

        est = np.array(
            [[[0.1, 0.02], [0.07, 0.25]], [[0.05, 0.1], [0.03, 0.025]]]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_eigen(self):
        err = self.m.compute_error(error_type="eigen", floor=False)

        est = np.array(
            [
                [[0.70890761, 0.70890761], [0.70890761, 0.70890761]],
                [[1.4186272, 1.4186272], [1.4186272, 1.4186272]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_eigen_floor(self):
        err = self.m.compute_error(error_type="eigen", floor=True)

        est = np.array(
            [
                [[0.70890761, 0.70890761], [0.70890761, 0.70890761]],
                [[1.4186272, 1.4186272], [1.4186272, 1.4186272]],
            ]
        )
        self.assertTrue(np.allclose(err, est))


# =============================================================================
# Run
# =============================================================================
if __name__ in "__main__":
    unittest.main()
