import os

from tests import TEST_MTPY_ROOT
import numpy as np
from mtpy.core.mt import MT

def test_rotate():
    alpha = 30
    alpharad = np.deg2rad(alpha)
    rotation_matrix = np.matrix([[np.cos(alpharad),np.sin(alpharad)],
                                 [-np.sin(alpharad),np.cos(alpharad)]])
    
    #mt_obj = MT(r'C:\mtpywin\mtpy\data\BBMT\EGC020A_pho.edi') # absolute path is user-dependent
    path2edi = os.path.normpath(os.path.join(TEST_MTPY_ROOT, 'data/BBMT/EGC020A_pho.edi'))
    mt_obj = MT(path2edi)

    z_array = mt_obj.Z.z.copy()
    
    ztest = z_array[0]
    
    ztest_rot = np.ma.dot(np.ma.dot(rotation_matrix,ztest),rotation_matrix.T)
    
    mt_obj.Z.rotate(alpha)
    
    assert np.all(mt_obj.Z.z[0] - ztest_rot < 1e-10)

    
if __name__ == "__main__":
    test_rotate()
