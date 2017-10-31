# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:19:35 2017

@author: u64125
"""

from mtpy.core.mt import MT
import mtpy.analysis.geometry as mtg
import numpy as np

mtObj = MT(r'C:\Git\mtpy\examples\data\edi_files\pb42c.edi')

strike_angle_result = np.array([[          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [          np.nan,           np.nan],
       [  38.45662316,  128.45662316],
       [  28.61883115,  118.61883115],
       [  14.45341494,  104.45341494],
       [   8.43320651,   98.43320651],
       [   4.94952784,   94.94952784],
       [   2.09090369,   92.09090369],
       [   1.39146887,   91.39146887],
       [   0.39905337,   90.39905337],
       [  -5.49553673,   84.50446327],
       [  -6.28846049,   83.71153951],
       [  -7.31641788,   82.68358212],
       [ -10.45341947,   79.54658053],
       [  -7.07075086,   82.92924914],
       [  -7.5429295 ,   82.4570705 ],
       [  -6.06405688,   83.93594312],
       [  -3.54915951,   86.45084049],
       [  -3.12596637,   86.87403363],
       [  -0.47404093,   89.52595907],
       [   2.74343665,   92.74343665],
       [   4.78078759,   94.78078759],
       [   7.71125988,   97.71125988],
       [  11.0123521 ,  101.0123521 ],
       [  13.81639678,  103.81639678],
       [  13.60497071,  103.60497071],
       [  15.87672806,  105.87672806]])
       
strike_angle = mtg.strike_angle(z_object = mtObj.Z)
   
assert np.all(np.abs(strike_angle[np.isfinite(strike_angle)] - \
                     strike_angle_result[np.isfinite(strike_angle_result)] < 1e-8))      
      