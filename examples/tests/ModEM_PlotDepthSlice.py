# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot Depth Slice 

"""
import os
import matplotlib.pyplot as plt

os.chdir(r'C:\Git\mtpy')
from mtpy.modeling.modem import PlotDepthSlice

# directory where files are located
wd = r'C:\Git\mtpy\examples\model_files\ModEM'

# directory to save to
save_path = r'C:\Git\mtpy\examples\plots\ModEM'

# file stem for inversion result
filestem = 'Modular_MPI_NLCG_004'

# period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
period_index = 0

# plot map
dsmap = PlotDepthSlice(model_fn = os.path.join(wd,filestem+'.rho'),
                        data_fn = os.path.join(wd,filestem+'dat'),
                        depth_index=30,
                        save_plots='n'
                        )
plt.savefig(os.path.join(save_path,'DepthSlice.png'))