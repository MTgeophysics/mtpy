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


class TestErrorEstimationZ(unittest.TestCase):
    def setUp(self):

        self.data = np.array(
            [
                [[0.1 - 0.1j, 10 + 10j], [-20 - 20j, -0.3 + 0.3j]],
                [[0.2 - 0.2j, 5 + 5j], [-40 - 40j, -0.6 + 0.6j]],
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
        est = np.array(
            [
                [[0.00707107, 0.70710678], [1.41421356, 0.0212132]],
                [[0.01414214, 0.35355339], [2.82842712, 0.04242641]],
            ]
        )

        self.assertTrue(np.allclose(err, est))

    def test_compute_percent_error_floor(self):
        err = self.m.compute_error(
            error_type="percent", error_value=5, floor=True
        )
        est = np.array(
            [
                [[0.1, 0.70710678], [1.41421356, 0.25]],
                [[0.05, 0.35355339], [2.82842712, 0.04242641]],
            ]
        )

        self.assertTrue(np.allclose(err, est))

    def test_compute_arithmetic_mean_error(self):
        err = self.m.compute_error(error_type="arithmetic_mean", floor=False)

        est = np.array(
            [
                [[1.06066017, 1.06066017], [1.06066017, 1.06066017]],
                [[1.59099026, 1.59099026], [1.59099026, 1.59099026]],
            ]
        )

        self.assertTrue(np.allclose(err, est))

    def test_compute_arithmetic_mean_error_floor(self):
        err = self.m.compute_error(error_type="arithmetic_mean", floor=True)

        est = np.array(
            [
                [[1.06066017, 1.06066017], [1.06066017, 1.06066017]],
                [[1.59099026, 1.59099026], [1.59099026, 1.59099026]],
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
                [[1.0001, 1.0001], [1.0001, 1.0001]],
                [[1.00039992, 1.00039992], [1.00039992, 1.00039992]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_eigen_floor(self):
        err = self.m.compute_error(error_type="eigen", floor=True)

        est = np.array(
            [
                [[1.0001, 1.0001], [1.0001, 1.0001]],
                [[1.00039992, 1.00039992], [1.00039992, 1.00039992]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_geometric_mean(self):
        err = self.m.compute_error(error_type="geometric_mean", floor=False)

        est = np.array([[[1.0, 1.0], [1.0, 1.0]], [[1.0, 1.0], [1.0, 1.0]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_geometric_mean_floor(self):
        err = self.m.compute_error(error_type="geometric_mean", floor=True)

        est = np.array([[[1.0, 1.0], [1.0, 1.0]], [[1.0, 1.0], [1.0, 1.0]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_by_row(self):
        err = self.m.compute_error(error_type="row", floor=False)

        est = np.array(
            [
                [[0.70710678, 0.70710678], [1.41421356, 1.41421356]],
                [[0.35355339, 0.35355339], [2.82842712, 2.82842712]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_by_row_floor(self):
        err = self.m.compute_error(error_type="row", floor=True)

        est = np.array(
            [
                [[0.70710678, 0.70710678], [1.41421356, 1.41421356]],
                [[0.35355339, 0.35355339], [2.82842712, 2.82842712]],
            ]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_absolute(self):
        err = self.m.compute_error(error_type="absolute", floor=False)

        est = np.array(
            [[[0.05, 0.05], [0.05, 0.05]], [[0.05, 0.05], [0.05, 0.05]]]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_absolute_floor(self):
        err = self.m.compute_error(error_type="absolute", floor=True)

        est = np.array(
            [[[0.1, 0.05], [0.07, 0.25]], [[0.05, 0.1], [0.05, 0.05]]]
        )
        self.assertTrue(np.allclose(err, est))


class TestErrorEstimationT(unittest.TestCase):
    def setUp(self):

        self.data = np.array(
            [[[0.25 - 0.2j, 0.25 + 0.2j]], [[0.15 - 0.3j, 0.75 + 0.8j]]]
        )
        self.data_error = np.array([[[0.03, 0.01]], [[0.05, 0.09]]])
        self.m = ModelErrors(self.data, self.data_error, mode="tipper")

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
        est = np.array(
            [[[0.01600781, 0.01600781]], [[0.01677051, 0.05482928]]]
        )

        self.assertTrue(np.allclose(err, est))

    def test_compute_percent_error_floor(self):
        err = self.m.compute_error(
            error_type="percent", error_value=5, floor=True
        )
        est = np.array([[[0.03, 0.01600781]], [[0.05, 0.09]]])

        self.assertTrue(np.allclose(err, est))

    def test_compute_arithmetic_mean_error(self):
        err = self.m.compute_error(error_type="arithmetic_mean", floor=False)

        est = np.array([[[0.01600781, 0.01600781]], [[0.0357999, 0.0357999]]])

        self.assertTrue(np.allclose(err, est))

    def test_compute_arithmetic_mean_error_floor(self):
        err = self.m.compute_error(error_type="arithmetic_mean", floor=True)

        est = np.array([[[0.03, 0.01600781]], [[0.05, 0.09]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_median(self):
        err = self.m.compute_error(error_type="median", floor=False)

        est = np.array([[[0.0125, 0.0125]], [[0.02573908, 0.02573908]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_median_floor(self):
        err = self.m.compute_error(error_type="median", floor=True)

        est = np.array([[[0.03, 0.0125]], [[0.05, 0.09]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_eigen(self):
        self.assertRaises(IndexError, self.m.compute_error, error_type="eigen")

    def test_compute_eigen_floor(self):
        self.assertRaises(
            IndexError, self.m.compute_error, error_type="eigen", floor=True
        )

    def test_compute_geometric_mean(self):
        err = self.m.compute_error(error_type="geometric_mean", floor=False)

        est = np.array(
            [[[0.01600781, 0.01600781]], [[0.03032351, 0.03032351]]]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_geometric_mean_floor(self):
        err = self.m.compute_error(error_type="geometric_mean", floor=True)

        est = np.array([[[0.03, 0.01600781]], [[0.05, 0.09]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_by_row(self):
        err = self.m.compute_error(error_type="row", floor=False)

        est = np.array(
            [[[0.01600781, 0.01600781]], [[0.01677051, 0.05482928]]]
        )
        self.assertTrue(np.allclose(err, est))

    def test_compute_by_row_floor(self):
        err = self.m.compute_error(error_type="row", floor=True)

        est = np.array([[[0.03, 0.01600781]], [[0.05, 0.09]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_absolute(self):
        err = self.m.compute_error(error_type="absolute", floor=False)

        est = np.array([[[0.05, 0.05]], [[0.05, 0.05]]])
        self.assertTrue(np.allclose(err, est))

    def test_compute_absolute_floor(self):
        err = self.m.compute_error(error_type="absolute", floor=True)

        est = np.array([[[0.05, 0.05]], [[0.05, 0.09]]])
        self.assertTrue(np.allclose(err, est))


# =============================================================================
# Run
# =============================================================================
if __name__ in "__main__":
    unittest.main()
