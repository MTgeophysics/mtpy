# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 13:46:49 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest
import numpy as np
import xarray as xr

from mtpy.core.transfer_function.base import TFBase

# =============================================================================


class TestTFBaseTFInput(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tf = TFBase(tf=np.array([[[0]], [[1]]]))
        self.expected_shape = (2, 1, 1)
        self.expected = {
            "transfer_function": {"dtype": complex},
            "transfer_function_error": {"dtype": float},
            "transfer_function_model_error": {"dtype": float},
        }

    def test_shape(self):
        for key, v_dict in self.expected.items():
            tf = getattr(self.tf._dataset, key)
            with self.subTest(f"{key} shape"):
                self.assertEqual(tf.shape, self.expected_shape)

            with self.subTest(f"{key} dtype"):
                self.assertEqual(tf.dtype, v_dict["dtype"])

    def test_zeros(self):
        with self.subTest("error"):
            self.assertEqual(
                (self.tf._dataset.transfer_function_error.values == 0).all(),
                True,
            )

        with self.subTest("model error"):
            self.assertEqual(
                (
                    self.tf._dataset.transfer_function_model_error.values == 0
                ).all(),
                True,
            )

    def test_frequency(self):
        self.assertEqual(
            (self.tf.frequency == 1.0 / np.arange(1, 3, 1)).all(), True
        )

    def test_period(self):
        self.assertEqual((self.tf.period == np.arange(1, 3, 1)).all(), True)


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
