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
            [[0.1 - 0.1j, 10 + 10j], [-0.3 + 0.3j, -20 - 20j]]
        )
        self.data_error = np.array([[0.1, 0.02], [0.07, 0.25]])
        self.m = ModelErrors(self.data, self.data_error)

    def test_bad_shape(self):
        def set_data(data):
            self.m.data = data

        self.assertRaises(ValueError, set_data, np.random.rand(3, 3, 3))

    def test_bad_shape_error(self):
        def set_measurement_error(error):
            self.m.measurement_error = error

        self.assertRaises(
            ValueError, set_measurement_error, np.random.rand(3, 3, 3)
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

    # def test_compute_percent_error(self):
    #     self.m.error_value = .


# =============================================================================
# Run
# =============================================================================
if __name__ in "__main__":
    unittest.main()
