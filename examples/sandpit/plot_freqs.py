"""
Description:
    what does this script module do? How to do it.

Author: fei.zhang@ga.gov.au

Date:
"""

__author__ = 'u25656'

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

freqs = [ 1.000000e+00,  5.623413e-01, 3.162278e-01 , 1.778279e-01 , 1.000000e-01 ,5.623413e-02,
 3.162278e-02 , 1.778279e-02 , 1.000000e-02 , 5.623413e-03 , 3.162278e-03 , 1.778279e-03,
 1.000000e-03,  5.623413e-04 , 3.162278e-04 , 1.778279e-04 , 1.000000e-04 ]

periods= [1/f for f in freqs]
plt.plot(freqs, '^', markersize='10')
plt.show()

plt.plot(periods, 'o', markersize='10')
plt.show()
