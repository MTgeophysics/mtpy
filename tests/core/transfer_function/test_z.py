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

    def test_z(self):
        self.assertEqual(None, self.z.z)

    def test_z_error(self):
        self.assertEqual(None, self.z.z_error)

    def test_z_model_error(self):
        self.assertEqual(None, self.z.z_model_error)

    def test_resistivity(self):
        self.assertEqual(None, self.z.resistivity)

    def test_phase(self):
        self.assertEqual(None, self.z.phase)

    def test_resistivity_error(self):
        self.assertEqual(None, self.z.resistivity_error)

    def test_phase_error(self):
        self.assertEqual(None, self.z.phase_error)

    def test_resistivity_model_error(self):
        self.assertEqual(None, self.z.resistivity_model_error)

    def test_phase_model_error(self):
        self.assertEqual(None, self.z.phase_model_error)


class TestZSetResPhase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.z = Z()
        self.resistivity = np.array([[[5.0, 100.0], [100.0, 5.0]]])
        self.phase = np.array([[[90.0, 45.0], [-135.0, -90.0]]])
        self.resistivity_error = np.array([[1, 5], [5, 1]])
        self.phase_error = np.array([[0.5, 1], [1, 0.5]])
        self.resistivity_model_error = np.array([[2, 10], [10, 2]])
        self.phase_model_error = np.array([[0.5, 1], [1, 0.5]])
        self.z.set_resistivity_phase(
            self.resistivity,
            self.phase,
            np.array([1]),
            res_error=self.resistivity_error,
            phase_error=self.phase_error,
            res_model_error=self.resistivity_model_error,
            phase_model_error=self.phase_model_error,
        )

    def test_is_empty(self):
        self.assertFalse(self.z._is_empty())

    def test_has_tf(self):
        self.assertTrue(self.z._has_tf())

    def test_has_tf_error(self):
        self.assertTrue(self.z._has_tf_error())

    def test_has_tf_model_error(self):
        self.assertTrue(self.z._has_tf_model_error())

    def test_resistivity(self):
        self.assertTrue(np.isclose(self.z.resistivity, self.resistivity).all())

    def test_phase(self):
        self.assertTrue(np.isclose(self.z.phase, self.phase).all())

    def test_resistivity_error(self):
        self.assertTrue(
            np.isclose(self.resistivity_error, self.z.resistivity_error).all()
        )

    def test_phase_error(self):
        self.assertTrue(
            np.isclose(
                self.z.phase_error,
                np.array(
                    [[[11.30993247, 2.86240523], [2.86240523, 11.30993247]]]
                ),
            ).all()
        )

    def test_resistivity_model_error(self):
        self.assertTrue(
            np.isclose(
                self.resistivity_model_error, self.z.resistivity_model_error
            ).all()
        )

    def test_phase_model_error(self):
        self.assertTrue(
            np.isclose(
                self.z.phase_model_error,
                np.array(
                    [[[21.80140949, 5.71059314], [5.71059314, 21.80140949]]]
                ),
            ).all()
        )


class TestRemoveStaticShift(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.z = Z()
        self.z.z = np.array(
            [[[0.1 - 0.1j, 10.0 + 10.0j], [-10.0 - 10.0j, -0.1 + 0.1j]]]
        )

        self.ss_x = 0.5
        self.ss_y = 1.5

        self.ss_z = np.array(
            [
                [
                    [0.14142136 - 0.14142136j, 14.14213562 + 14.14213562j],
                    [-8.16496581 - 8.16496581j, -0.08164966 + 0.08164966j],
                ]
            ]
        )

        self.new_z = self.z.remove_ss(self.ss_x, self.ss_y, inplace=False)

    def test_remove_ss(self):
        self.assertTrue(np.allclose(self.new_z.z, self.ss_z))

    def test_ss_factors(self):
        self.assertTrue(
            np.allclose(
                (self.z.z / self.new_z.z) ** 2,
                np.array(
                    [[[0.5 + 0.0j, 0.5 + 0.0j], [1.5 - 0.0j, 1.5 - 0.0j]]]
                ),
            )
        )

    def test_set_factor_fail_single(self):
        self.assertRaises(ValueError, self.z.remove_ss, "k")

    def test_set_factor_fail_too_many(self):
        self.assertRaises(ValueError, self.z.remove_ss, [1, 2, 3])


class TestRemoveDistortion(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        z = np.array(
            [[[0.1 - 0.1j, 10.0 + 10.0j], [-10.0 - 10.0j, -0.1 + 0.1j]]]
        )
        self.z = Z()
        self.z.z = z

        self.distortion = np.matrix(np.real([[0.5, 0.1], [0.2, 0.6]]))

        self.dz = Z()
        self.dz.z = np.array(np.dot(self.distortion, z)).reshape((1, 2, 2))

    def test_remove_distortion(self):
        d, new_z = self.dz.remove_distortion(self.distortion)

        with self.subTest("distortion matrix"):
            self.assertTrue(np.allclose(self.distortion, d))

        with self.subTest(("z")):
            self.assertTrue(np.allclose(new_z.z, self.z.z))

    def test_fail_bad_input_shape_too_many(self):
        self.assertRaises(
            ValueError, self.z.remove_distortion, np.random.rand(4, 3, 3, 3)
        )

    def test_fail_bad_input_shape_not_z_shape(self):
        self.assertRaises(
            ValueError, self.z.remove_distortion, np.random.rand(4, 3, 3)
        )

    def test_fail_bad_input_singular_matrix(self):
        self.assertRaises(
            ValueError, self.z.remove_distortion, np.matrix([[0, 1], [0, 0]])
        )


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
