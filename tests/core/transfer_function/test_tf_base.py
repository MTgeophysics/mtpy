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
import scipy.signal as spi
from mtpy.core.transfer_function.base import TFBase
from mtpy.utils.calculator import (
    rotate_matrix_with_errors,
    rotate_vector_with_errors,
)

# =============================================================================


class TestTFBaseTFInput(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tf = TFBase(tf=np.array([[[0, 1], [1, 0]], [[1, 0], [0, 1]]]))
        self.expected_shape = (2, 2, 2)
        self.expected = {
            "transfer_function": {"dtype": complex, "empty": False},
            "transfer_function_error": {"dtype": float, "empty": True},
            "transfer_function_model_error": {"dtype": float, "empty": True},
        }

    def test_shape_zeros_dtype(self):
        for key, v_dict in self.expected.items():
            tf = getattr(self.tf._dataset, key)
            with self.subTest(f"{key} shape"):
                self.assertEqual(tf.shape, self.expected_shape)

            with self.subTest(f"{key} dtype"):
                self.assertEqual(tf.dtype, v_dict["dtype"])

            with self.subTest(f"{key} empty"):
                self.assertEqual((tf.values == 0).all(), v_dict["empty"])

    def test_frequency(self):
        self.assertEqual(
            (self.tf.frequency == 1.0 / np.arange(1, 3, 1)).all(), True
        )

    def test_period(self):
        self.assertEqual((self.tf.period == np.arange(1, 3, 1)).all(), True)


class TestTFBaseTFErrorInput(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tf = TFBase(
            tf_error=np.array([[[0, 1], [1, 0]], [[1, 0], [0, 1]]])
        )
        self.expected_shape = (2, 2, 2)
        self.expected = {
            "transfer_function": {"dtype": complex, "empty": True},
            "transfer_function_error": {"dtype": float, "empty": False},
            "transfer_function_model_error": {"dtype": float, "empty": True},
        }

        def test_shape_zeros_dtype(self):
            for key, v_dict in self.expected.items():
                tf = getattr(self.tf._dataset, key)
                with self.subTest(f"{key} shape"):
                    self.assertEqual(tf.shape, self.expected_shape)

                with self.subTest(f"{key} dtype"):
                    self.assertEqual(tf.dtype, v_dict["dtype"])

                with self.subTest(f"{key} empty"):
                    self.assertEqual((tf.values == 0).all(), v_dict["empty"])

    def test_frequency(self):
        self.assertEqual(
            (self.tf.frequency == 1.0 / np.arange(1, 3, 1)).all(), True
        )

    def test_period(self):
        self.assertEqual((self.tf.period == np.arange(1, 3, 1)).all(), True)


class TestTFBaseTFModelErrorInput(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tf = TFBase(
            tf_model_error=np.array([[[0, 1], [1, 0]], [[1, 0], [0, 1]]])
        )
        self.expected_shape = (2, 1, 1)
        self.expected = {
            "transfer_function": {"dtype": complex, "empty": True},
            "transfer_function_error": {"dtype": float, "empty": True},
            "transfer_function_model_error": {"dtype": float, "empty": False},
        }

        def test_shape_zeros_dtype(self):
            for key, v_dict in self.expected.items():
                tf = getattr(self.tf._dataset, key)
                with self.subTest(f"{key} shape"):
                    self.assertEqual(tf.shape, self.expected_shape)

                with self.subTest(f"{key} dtype"):
                    self.assertEqual(tf.dtype, v_dict["dtype"])

                with self.subTest(f"{key} empty"):
                    self.assertEqual((tf.values == 0).all(), v_dict["empty"])

    def test_frequency(self):
        self.assertEqual(
            (self.tf.frequency == 1.0 / np.arange(1, 3, 1)).all(), True
        )

    def test_period(self):
        self.assertEqual((self.tf.period == np.arange(1, 3, 1)).all(), True)


class TestTFBaseFrequencyInput(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tf = TFBase(frequency=[1, 2, 3])
        self.expected_shape = (3, 2, 2)
        self.expected = {
            "transfer_function": {"dtype": complex, "empty": True},
            "transfer_function_error": {"dtype": float, "empty": True},
            "transfer_function_model_error": {"dtype": float, "empty": True},
        }

    def test_set_frequency(self):
        self.tf.frequency = np.logspace(-1, 1, 3)
        with self.subTest("freq"):
            self.assertEqual(
                (self.tf.frequency == np.logspace(-1, 1, 3)).all(), True
            )
        with self.subTest("period"):
            self.assertEqual(
                (self.tf.period == 1.0 / np.logspace(-1, 1, 3)).all(), True
            )

    def test_set_period(self):
        self.tf.period = 1.0 / np.logspace(-1, 1, 3)
        with self.subTest("freq"):
            self.assertEqual(
                (self.tf.frequency == np.logspace(-1, 1, 3)).all(), True
            )
        with self.subTest("period"):
            self.assertEqual(
                (self.tf.period == 1.0 / np.logspace(-1, 1, 3)).all(), True
            )


class TestTFBaseValidators(unittest.TestCase):
    def setUp(self):
        self.tf = TFBase()

    def test_validate_array_input_float(self):
        self.assertEqual(
            (
                np.zeros((1, 2, 2), dtype=float)
                == self.tf._validate_array_input(
                    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]], float
                )
            ).all(),
            True,
        )

    def test_validate_array_input_complex(self):
        self.assertEqual(
            (
                np.zeros((1, 2, 2), dtype=complex)
                == self.tf._validate_array_input(
                    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]], complex
                )
            ).all(),
            True,
        )

    def test_validate_array_input_int(self):
        self.assertEqual(
            (
                np.zeros((1, 2, 2), dtype=float)
                == self.tf._validate_array_input(
                    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]], float
                )
            ).all(),
            True,
        )

    def test_validate_frequency_shape(self):
        self.assertEqual(self.tf._validate_frequency([1], 10).size, 10)

    def test_is_empty(self):
        self.assertEqual(self.tf._is_empty(), True)

    def test_has_tf(self):
        self.assertEqual(self.tf._has_tf(), False)

    def test_has_tf_error(self):
        self.assertEqual(self.tf._has_tf_error(), False)

    def test_has_tf_model_error(self):
        self.assertEqual(self.tf._has_tf_model_error(), False)


class TestTFRotation(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tf = TFBase(
            tf=np.ones((3, 2, 2)),
            tf_error=np.ones((3, 2, 2)) * 0.25,
            tf_model_error=np.ones((3, 2, 2)) * 0.5,
        )

        self.rot_tf = self.tf.rotate(30)

        self.true_rot_tf = np.zeros((3, 2, 2), dtype=complex)
        self.true_rot_tf_error = np.zeros((3, 2, 2), dtype=float)
        self.true_rot_tf_model_error = np.zeros((3, 2, 2), dtype=float)
        for ii, angle in enumerate([30, 30, 30]):
            (
                self.true_rot_tf[ii],
                self.true_rot_tf_error[ii],
            ) = rotate_matrix_with_errors(
                np.ones((2, 2), dtype=complex), 30, np.ones((2, 2)) * 0.25
            )
            (_, self.true_rot_tf_model_error[ii],) = rotate_matrix_with_errors(
                np.ones((2, 2), dtype=complex), 30, np.ones((2, 2)) * 0.5
            )

    def test_tf(self):
        self.assertEqual(
            (
                self.tf._dataset.transfer_function.values
                == np.ones((3, 2, 2), dtype=complex)
            ).all(),
            True,
        )

    def test_rot_tf(self):
        self.assertEqual(
            np.isclose(
                self.rot_tf.transfer_function.values, self.true_rot_tf
            ).all(),
            True,
        )

    def test_tf_error(self):
        self.assertEqual(
            (
                self.tf._dataset.transfer_function_error.values
                == np.ones((3, 2, 2), dtype=float) * 0.25
            ).all(),
            True,
        )

    def test_rot_tf_error(self):
        self.assertEqual(
            np.isclose(
                self.rot_tf.transfer_function_error.values,
                self.true_rot_tf_error,
            ).all(),
            True,
        )

    def test_tf_model_error(self):
        self.assertEqual(
            (
                self.tf._dataset.transfer_function_model_error.values
                == np.ones((3, 2, 2), dtype=float) * 0.5
            ).all(),
            True,
        )

    def test_rot_tf_model_error(self):
        self.assertEqual(
            np.isclose(
                self.rot_tf.transfer_function_model_error.values,
                self.true_rot_tf_model_error,
            ).all(),
            True,
        )


class TestTFInterpolation(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.period = np.logspace(-3, 3, 6)
        self.t = np.linspace(0, 24, 24) * 0.1
        self.tf = np.array(
            [
                np.cos(pp * np.pi * 2 * 10 * self.t)
                + 1j * np.sin(pp * np.pi * 2 * 10 * self.t)
                for pp in self.period
            ]
        ).sum(axis=0)

        self.tf = self.tf.reshape((6, 2, 2))
        self.tf_error = np.abs(self.tf) * 0.05
        self.tf_model_error = np.abs(self.tf) * 0.10

        self.tf_base = TFBase(
            tf=self.tf,
            tf_error=self.tf_error,
            tf_model_error=self.tf_model_error,
            frequency=1.0 / self.period,
        )

        self.new_periods = np.logspace(-3, 3, 12)

    def interpolate(self, interp_type, bounds_error=False):
        zmap = {0: "x", 1: "y"}
        interp_dict = {}

        interp_tf = np.zeros((self.new_periods.size, 2, 2), dtype=complex)
        interp_tf_error = np.zeros((self.new_periods.size, 2, 2), dtype=float)
        interp_tf_model_error = np.zeros(
            (self.new_periods.size, 2, 2), dtype=float
        )
        for ii in range(2):
            for jj in range(2):
                comp = f"{zmap[ii]}{zmap[jj]}"
                interp_dict[comp] = {}
                # need to look out for zeros in the impedance
                # get the indicies of non-zero components
                nz_index = np.nonzero(self.Z.z[:, ii, jj])

                if len(nz_index[0]) == 0:
                    continue
                # get the non-zero components
                tf_real = self.tf[nz_index, ii, jj].real
                tf_imag = self.tf[nz_index, ii, jj].imag
                tf_err = self.tf_error[nz_index, ii, jj]
                tf_model_err = self.tf_error[nz_index, ii, jj]

                # get the frequencies of non-zero components
                f = 1.0 / self.period

                # create a function that does 1d interpolation
                tr = spi.interp1d(
                    f,
                    tf_real,
                    kind=interp_type,
                    bounds_error=bounds_error,
                )
                ti = spi.interp1d(
                    f,
                    tf_imag,
                    kind=interp_type,
                    bounds_error=bounds_error,
                )
                te = spi.interp1d(
                    f,
                    tf_err,
                    kind=interp_type,
                    bounds_error=bounds_error,
                )
                tme = spi.interp1d(
                    f,
                    tf_model_err,
                    kind=interp_type,
                    bounds_error=bounds_error,
                )

                interp_tf[:, ii, jj] = tr(self.new_periods) + 1j * ti(
                    self.new_periods
                )
                interp_tf_error[:, ii, jj] = te(self.new_periods)
                interp_tf_model_error[:, ii, jj] = tme(self.new_periods)

        interp_ds = TFBase(
            tf=interp_tf,
            tf_error=interp_tf_error,
            tf_model_error=interp_tf_model_error,
            period=self.new_periods,
        )

        return interp_ds

    def test_nearest(self):
        interp_ds = self.interpolate("nearest")
        interp_tf = self.tf_base.interpolate(self.new_periods, method="nearest")

        for key in [
            "transfer_function",
            "transfer_function_error",
            "transfer_function_model_error",
        ]:

            with self.subTest(key):
                self.assertEqual(
                    np.isclose(
                        interp_ds._dataset[key].values,
                        interp_tf._dataset[key].values,
                    ),
                    True,
                )


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
