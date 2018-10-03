from unittest import TestCase
import numpy as np
import pytest

from mtpy.utils.calculator import get_period_list, make_log_increasing_array,\
                                  z_error2r_phi_error, nearest_index


class TestCalculator(TestCase):
    def setUp(self):
        
        self.z1_layer = 80
        self.target_depth = 400e3
        self.n_layers = 120
        self.increment_factor = 0.999
        self.z = np.array([[[ -6.4-1.320e+01j,  45.2+9.920e+01j],
                            [-42.5-1.014e+02j,  10.4+1.930e+01j]],

                           [[ -1. -2.100e+00j,   8.7+1.590e+01j],
                            [ -7.5-1.520e+01j,   1.7+3.100e+00j]],
                    
                           [[ -1.3-1.000e-01j,   6.3+2.000e+00j],
                            [ -6.3-1.600e+00j,   1.6+4.000e-01j]]])
        self.z_err = np.array([[[ 1.5, 10.9],
                                [11. ,  2.2]],
                        
                               [[ 0.2,  1.8],
                                [ 1.7,  0.4]],
                        
                               [[ 0.1,  0.7],
                                [ 0.6,  2.0]]])
        self.freq = np.array([100.,10.,1.])
        
        return
        
    
    def test_nearest_index(self):
        
        freqfind = 8.
        self.assertTrue(nearest_index(freqfind,self.freq)==1)
        
        freqfind2 = 1.2
        self.assertTrue(nearest_index(freqfind2,self.freq)==2)
        
        freqfind3 = 1000.
        self.assertTrue(nearest_index(freqfind3,self.freq)==0)


    def test_get_period_list(self):
        
        array1 = np.array([1.77827941e-02, 3.16227766e-02, 5.62341325e-02, 1.00000000e-01,
                           1.77827941e-01, 3.16227766e-01, 5.62341325e-01, 1.00000000e+00,
                           1.77827941e+00, 3.16227766e+00, 5.62341325e+00, 1.00000000e+01,
                           1.77827941e+01, 3.16227766e+01, 5.62341325e+01, 1.00000000e+02,
                           1.77827941e+02, 3.16227766e+02, 5.62341325e+02])
        array2 = np.array([3.16227766e-02, 5.62341325e-02, 1.00000000e-01,
                           1.77827941e-01, 3.16227766e-01, 5.62341325e-01, 1.00000000e+00,
                           1.77827941e+00, 3.16227766e+00, 5.62341325e+00, 1.00000000e+01,
                           1.77827941e+01, 3.16227766e+01, 5.62341325e+01, 1.00000000e+02,
                           1.77827941e+02, 3.16227766e+02])
        
        array1test = get_period_list(0.02,400,4,include_outside_range=True)
        array2test = get_period_list(0.02,400,4,include_outside_range=False)
        
        self.assertTrue(all(np.abs(array1test-array1)/array1 < 1e-8))
        self.assertTrue(all(np.abs(array2test-array2)/array2 < 1e-8))
        
        
        # test with ends of input range on an exact decade
        array1 = np.array([  0.1       ,   0.21544347,   0.46415888,   1.        ,
                           2.15443469,   4.64158883,  10.        ,  21.5443469 ,
                           46.41588834, 100.        ])
        
        array1test = get_period_list(0.1,100,3,include_outside_range=True)
        array2test = get_period_list(0.1,100,3,include_outside_range=False)
        
        self.assertTrue(all(np.abs(array1test-array1)/array1 < 1e-8))
        self.assertTrue(all(np.abs(array2test-array1)/array1 < 1e-8))
        
    
    def test_make_log_increasing_array(self):
        
        array1 = make_log_increasing_array(self.z1_layer,
                                           self.target_depth,
                                           self.n_layers,
                                           increment_factor = self.increment_factor)
        
        self.assertTrue(np.abs(array1.sum()/self.target_depth - 1.)\
                        < 1. - self.increment_factor)
        
    
    def test_z_error2r_phi_error(self):
        # relative error in resistivity is 2 * relative error in z
        res_rel_err_test = 2. * self.z_err/np.abs(self.z)

        phase_err_test = np.rad2deg(np.arctan(self.z_err/np.abs(self.z)))
        phase_err_test[res_rel_err_test > 1.] = 90.
        
        # test providing an array
        res_rel_err, phase_err = z_error2r_phi_error(self.z.real,self.z.imag, self.z_err)
        
        self.assertTrue(np.all(np.abs(res_rel_err-res_rel_err_test)/res_rel_err_test < 1e-8))
        self.assertTrue(np.all(np.abs(phase_err-phase_err_test)/phase_err_test < 1e-8))
        
        # test providing a single value
        res_rel_err, phase_err = z_error2r_phi_error(self.z.real[0,0,1],self.z.imag[0,0,1], self.z_err[0,0,1])
        
        self.assertTrue(np.all(np.abs(res_rel_err-res_rel_err_test[0,0,1])/res_rel_err_test[0,0,1] < 1e-8))
        self.assertTrue(np.all(np.abs(phase_err-phase_err_test[0,0,1])/phase_err_test[0,0,1] < 1e-8))        
