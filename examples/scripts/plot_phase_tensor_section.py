# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period) 
"""

from mtpy.imaging.phase_tensor_pseudosection import PlotPhaseTensorPseudoSection
import os.path as op
import os
import matplotlib.pyplot as plt

# path to edis
edi_path = r'C:\mtpywin\mtpy\examples\data\edi_files_2'

# edi list
elst=[op.join(edi_path,edi) for edi in os.listdir(edi_path) if ((edi.endswith('.edi')))]# and edi.startswith('GB')

# create a plot object
plotObj = PlotPhaseTensorPseudoSection(fn_list = elst,
                                 linedir='ns', # 'ns' if the line is closer to north-south, 'ew' if line is closer to east-west
                                 stretch=(35,8), # determines (x,y) aspect ratio of plot
                                 station_id=(0,10), # indices for showing station names
                                 plot_tipper = 'yri', # plot tipper ('y') + 'ri' means real+imag
                                 font_size=5,
                                 lw=0.5
#                                 dpi=300
                                 )

# update some parameters
plotObj.ellipse_size = 2.5

## example to color by skew
#plotObj.ellipse_colorby = 'skew'
#plotObj.ellipse_cmap = 'mt_seg_bl2wh2rd'
#plotObj.ellipse_range = (-12,12,3)

# example to color by phimin
plotObj.ellipse_colorby = 'phimin'



plotObj.plot()