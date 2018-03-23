import numpy as np
from mtpy.core.mt import MT

def test_rotate():
    alpha = 30
    alpharad = np.deg2rad(alpha)
    rotation_matrix = np.matrix([[np.cos(alpharad),np.sin(alpharad)],
                                 [-np.sin(alpharad),np.cos(alpharad)]])
    
    mt_obj = MT(r'C:\mtpywin\mtpy\data\BBMT\EGC020A_pho.edi')
    
    z_array = mt_obj.Z.z.copy()
    
    ztest = z_array[0]
    
    ztest_rot = np.ma.dot(np.ma.dot(rotation_matrix,ztest),rotation_matrix.T)
    
    mt_obj.Z.rotate(alpha)
    
    assert np.all(mt_obj.Z.z[0] - ztest_rot < 1e-10)
    
    
    
if __name__ == "__main__":
    test_rotate()
