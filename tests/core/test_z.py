"""
Test mtpy.core.z

Updated: Peacock (11/2020)
    - Adding tests for inputs and coherence
 
"""
# =============================================================================
# 
# =============================================================================
import unittest
import numpy as np
from pathlib import Path

# import pytest
from tests import TEST_MTPY_ROOT
from mtpy.core.z import Z
from mtpy.core.mt import MT
from mtpy.utils.exceptions import MTpyError_Z

# =============================================================================
# Test only Z
# =============================================================================
class TestOnlyZ(unittest.TestCase):
    def setUp(self):
        self.z_obj = Z()
        
    def test_2x2_input(self):
        self.z_obj.z = np.random.random((2, 2)) + 1j * np.random.random((2, 2))
        
        self.assertEqual(self.z_obj.z.shape, (1, 2, 2))
        self.assertIn(self.z_obj.z.dtype, ["complex"])
        
    def test_fail_z_input(self):
        self.assertRaises(MTpyError_Z, self.z_obj.z, [[1], [2], [3]])
        
    def test_list_input(self):
        self.z_obj.z = [[8, 9], [10, 11]]
        
        self.assertEqual(self.z_obj.z.shape, (1, 2, 2))
        self.assertIn(self.z_obj.z.dtype, ["complex"])
        
    def test_freq_input(self):
        
        self.z_obj.freq = [1, 2, 3]
        
        self.assertIsInstance(self.z_obj.freq, np.ndarray)
        
    def test_z_err_input(self):
        self.z_obj.z_err = [[1, 2], [3, 4]]
        
        self.assertEqual(self.z_obj.z_err.shape, (1, 2, 2))
        self.assertIn(self.z_obj.z_err.dtype, ["float"])
        
    def test_fail_input_lengths(self):
        self.z_obj.z = np.random.random((5, 2, 2)) + 1j * np.random.random((5, 2, 2))
        
        self.assertRaises(MTpyError_Z, self.z_obj.freq, [1, 2])
        self.assertRaises(MTpyError_Z, self.z_obj.z_err, np.random.random((2, 2, 2)))
        
    def test_res_phase_calc(self):
        self.z_obj.z = np.random.random((5, 2, 2)) + 1j * np.random.random((5, 2, 2))
        self.z_obj.freq = np.logspace(0, 3, 5)
        self.z_obj.z_err = np.random.random((5, 2, 2))
        
        self.assertNotIsInstance(self.z_obj.res_xx, type(None))
        self.assertNotIsInstance(self.z_obj.phase_xx, type(None))
        self.assertNotIsInstance(self.z_obj.res_det, type(None))
        self.assertNotIsInstance(self.z_obj.res_err_xx, type(None))
        self.assertNotIsInstance(self.z_obj.phase_err_xx, type(None))
        self.assertNotIsInstance(self.z_obj.res_det_err, type(None))

# =============================================================================
# Test input of an edi file
# =============================================================================
class TestZ(unittest.TestCase):
    def setUp(self):
        edi_fn = Path(TEST_MTPY_ROOT, r"data/BBMT/EGC020A_pho.edi")
        self.MT = MT(edi_fn)
        self.rotation_angle = 30
        self.static_shift_x = 1.2
        self.static_shift_y = 1.5

    def test_rotate(self):

        alpharad = np.deg2rad(self.rotation_angle)
        rotation_matrix = np.matrix(
            [
                [np.cos(alpharad), np.sin(alpharad)],
                [-np.sin(alpharad), np.cos(alpharad)],
            ]
        )

        z_array = self.MT.Z.z.copy()

        ztest = z_array[0]

        ztest_rot = np.ma.dot(np.ma.dot(rotation_matrix, ztest), rotation_matrix.T)

        self.MT.Z.rotate(self.rotation_angle)

        self.assertTrue(np.all(self.MT.Z.z[0] - ztest_rot < 1e-10))

    def test_remove_ss(self):
        # calculate a corrected z array
        zcor = self.MT.Z.remove_ss(
            reduce_res_factor_x=self.static_shift_x,
            reduce_res_factor_y=self.static_shift_y,
        )[1]

        # make a Z object with the corrected array and calculate resistivity and phase
        Zcor = Z(z_array=zcor, freq=self.MT.Z.freq)
        Zcor.compute_resistivity_phase()

        # check that the resistivity factors are correct
        self.assertTrue(
            np.all(
                np.array((self.MT.Z.resistivity / Zcor.resistivity)[:, 0])
                - self.static_shift_x
                < 1e10
            )
        )
        self.assertTrue(
            np.all(
                np.array((self.MT.Z.resistivity / Zcor.resistivity)[:, 1])
                - self.static_shift_y
                < 1e10
            )
        )

    def test_compute_resistivity_phase(self):

        freq = np.logspace(2, -3, 6)
        z = np.array(
            [
                [
                    [-6.395642 - 1.316922e01j, 45.22294 + 9.921218e01j],
                    [-42.52678 - 1.013927e02j, 10.36037 + 1.929462e01j],
                ],
                [
                    [-1.028196 - 2.087766e00j, 8.6883 + 1.588160e01j],
                    [-7.535859 - 1.517449e01j, 1.724745 + 3.063922e00j],
                ],
                [
                    [-1.255376 - 7.618381e-02j, 6.325392 + 1.997068e00j],
                    [-6.281115 - 1.554746e00j, 1.635191 + 4.083652e-01j],
                ],
                [
                    [-1.285838 - 1.969697e-01j, 5.428554 + 1.687451e00j],
                    [-5.727512 - 1.275399e00j, 1.083837 + 4.377386e-01j],
                ],
                [
                    [-1.065314 - 4.287129e-01j, 2.326234 + 1.113632e00j],
                    [-2.71029 - 2.188917e00j, 0.6587484 + 7.972820e-02j],
                ],
                [
                    [-0.2464714 - 3.519637e-01j, 1.588412 + 5.252950e-01j],
                    [-0.6028374 - 8.850180e-01j, 0.5539547 + 1.885912e-01j],
                ],
            ]
        )
        zObj = Z(z_array=z, freq=freq)
        zObj.compute_resistivity_phase()

        # resistivity array computed from above z and freq arrays and verified
        # against outputs from WinGLink
        res_test = np.array(
            [
                [[4.286652e-01, 2.377634e01], [2.417802e01, 9.592394e-01]],
                [[1.083191e-01, 6.554236e00], [5.741087e00, 2.472473e-01]],
                [[3.163546e-01, 8.799773e00], [8.373929e00, 5.681224e-01]],
                [[3.384353e00, 6.463338e01], [6.886208e01, 2.732636e00]],
                [[2.637378e01, 1.330308e02], [2.427406e02, 8.806122e00]],
                [[3.692532e01, 5.597976e02], [2.293340e02, 6.848650e01]],
            ]
        )

        self.assertTrue(np.all(np.abs(zObj.resistivity / res_test - 1.0) < 1e-6))

# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()