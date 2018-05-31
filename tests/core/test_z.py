from unittest import TestCase
import numpy as np
#import pytest
from tests import TEST_MTPY_ROOT
from mtpy.core.z import Z
from mtpy.core.mt import MT
import os

class TestZ(TestCase):
    def setUp(self):
        self.edi_file = os.path.join(TEST_MTPY_ROOT, r'data/BBMT/EGC020A_pho.edi')
        self.MT = MT(self.edi_file)
        self.rotation_angle = 30
        self.static_shift_x = 1.2
        self.static_shift_y = 1.5


    def test_rotate(self):

        alpharad = np.deg2rad(self.rotation_angle)
        rotation_matrix = np.matrix([[np.cos(alpharad),np.sin(alpharad)],
                                     [-np.sin(alpharad),np.cos(alpharad)]])
        
        z_array = self.MT.Z.z.copy()
        
        ztest = z_array[0]
        
        ztest_rot = np.ma.dot(np.ma.dot(rotation_matrix,ztest),rotation_matrix.T)
        
        self.MT.Z.rotate(self.rotation_angle)
        
        self.assertTrue(np.all(self.MT.Z.z[0] - ztest_rot < 1e-10))


    def test_remove_ss(self):
        # calculate a corrected z array
        zcor = self.MT.Z.remove_ss(reduce_res_factor_x=self.static_shift_x,
                                   reduce_res_factor_y=self.static_shift_y)[1]
        
        # make a Z object with the corrected array and calculate resistivity and phase
        Zcor = Z(z_array = zcor, freq=self.MT.Z.freq)
        Zcor.compute_resistivity_phase()
        
        # check that the resistivity factors are correct
        self.assertTrue(np.all(np.array((self.MT.Z.resistivity/Zcor.resistivity)[:,0]) - self.static_shift_x < 1e10))
        self.assertTrue(np.all(np.array((self.MT.Z.resistivity/Zcor.resistivity)[:,1]) - self.static_shift_y < 1e10))