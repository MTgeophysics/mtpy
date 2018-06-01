from unittest import TestCase
import numpy as np
import pytest

from mtpy.utils.calculator import get_period_list, make_log_increasing_array


class TestCalculator(TestCase):
    def setUp(self):
        
        self.z1_layer = 80
        self.target_depth = 400e3
        self.n_layers = 120
        self.increment_factor = 0.999

        
        return
        


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
        
    
