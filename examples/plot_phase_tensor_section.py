# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period) 
"""

import mtpy.imaging.plotptpseudosection as ptp
import os.path as op
import os
import matplotlib.pyplot as plt

# path to edis
epath = r'C:\Git\mtpy\examples\data\edi_files\georgina'


elst=[op.join(epath,edi) for edi in os.listdir(epath) if edi.endswith('.edi')]


ptp.PlotPhaseTensorPseudoSection(fn_list = elst,
                                 tscale = 'period',
                                 ylim = (1e-1,1e3), # period range to plot
                                 #xlim = (0,10000),
                                 stretch=(50,16), # determines (x,y) aspect ratio of plot
                                 station_id=(0,10), # indices for showing station names
                                 ellipse_dict={'size':3},
                                 plot_tipper = 'yr',
                                 arrow_dict = {'size':3,'head_length':0.1,
                                               'head_width':0.1,'lw':0.5},# arrow parameters, adjust as necessary. lw = linewidth
                                 font_size=4,
                                 dpi=300)
