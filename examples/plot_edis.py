# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period) 
"""

import mtpy.imaging.plotresponse as mtpr
import mtpy.imaging.plotptpseudosection as ppts
import mtpy.core.edi as mtedi
import mtpy.analysis.geometry as mtg
import os.path as op
import os
import matplotlib.pyplot as plt
import numpy as np



epath = r'C:\Git\mtpy\examples\data\edi_files\georgina'

svdir = r'V:\2008\09GA_GA1_Georgina\Reprocessing_2016\plots\station_plots'

elst=[op.join(epath,edi) for edi in os.listdir(epath) if (edi.endswith('.edi'))]# and edi.startswith('GB')

dimensionality = np.zeros((len(elst),46))


for efile in elst:
#    eo = mtedi.Edi(filename=efile)
	pr = mtpr.PlotResponse(fn=efile,plot_num=2,res_limits=(1,10000),phase_limits=(0,90))
	plt.close()

