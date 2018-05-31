from tests import TEST_MTPY_ROOT
import numpy as np
from mtpy.utils.calculator import get_logspace_array

def test_get_logspace_array():
    
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
    
    array1test = get_logspace_array(0.02,400,4,include_outside_range=True)
    array2test = get_logspace_array(0.02,400,4,include_outside_range=False)
    
    assert all(np.abs(array1test-array1)/array1 < 1e-8)
    assert all(np.abs(array2test-array2)/array2 < 1e-8)
    
if __name__ == "__main__":
    test_get_logspace_array()
