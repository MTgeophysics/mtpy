# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:43:39 2021

@author: u64125
"""

from mtpy.imaging.plot_resphase_maps import PlotResPhaseMaps
import os

edipath = r'C:\mtpywin\mtpy\examples\data\ET_edi'


# frequency to plot
plot_freq = 1e-2


# gets edi file names as a list
elst = [os.path.join(edipath, f) for f in os.listdir(edipath) if f.endswith('.edi')]

RPMap = PlotResPhaseMaps(fn_list=elst)



# ### OR...  to use MT object list
# from mtpy.core.mt import MT
# mt_obj_list = [MT(fn) for fn in elst]


# RPMap = PlotResPhaseMaps(mt_list=mt_obj_list)




# make plot
RPMap.plot(plot_freq,
           'res',
           1,
           1000,
           extrapolation_buffer_degrees=0.2
           )