# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plots data/model misfit for a given period, at all sites, separately for each
of the impedance tensor modes + tipper

"""


import os
import os.path as op
import numpy as np

os.chdir(r'C:\mtpywin\mtpy')

from mtpy.modeling.modem import PlotRMSMaps
#from mtpy.imaging.plot_rms_map import PlotRMSMaps


wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'
savepath = r'C:\tmp'



filestem = op.join(wd,'Modular_MPI_NLCG_004')
resid_fn=op.join(wd,filestem + '.res')

probj = PlotRMSMaps(resid_fn,
                    period_index='all',
                    rms_cmap='jet', 
                    rms_max=5
                    )
probj.save_figure(save_path=savepath,
                  save_fig_dpi = 400 # change to your preferred figure resolution
                  )