# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 07:29:58 2013

@author: Alison Kirkby

plots phase tensor ellipses as a map for a given frequency
"""

import mtpy.imaging.mtplottools as mtpl
import os
import os.path as op
import matplotlib.pyplot as plt
import mtpy.imaging.plotptmaps as pptmaps


# directory containing edis
#edipath = r'C:\Git\mtpy\examples\data\edi_files\georgina'
edipath = r'/Softlab/Githubz/mtpy2/examples/data/edi_files/georgina'

# whether or not to save the figure to file
save = True

# full path to file to save to
#savepath = r'C:\Git\mtpy\examples\data\edi_files\georgina'
savepath = r'/tmp/georgina'

# frequency to plot
plot_freq = 0.1

# gets edi file names as a list
elst = [op.join(edipath,f) for f in os.listdir(edipath) if f.endswith('.edi')]

# parameters describing ellipses
ellipse_dict = {'size':0.02,'colorby':'phimin','range':(0,90,1),'cmap':'mt_bl2gr2rd'}
                     
# parameters describing the induction vector arrows
arrow_dict = {'size':0.2,
              'lw':0.01,
              'head_width':0.002,
              'head_length':0.002,
              'direction':0,
              'threshold':0.8,
              'direction':0}
              
              
# parameters describing the arrow legend (should be self explanatory)
arrow_legend_dict = {'position':'upper right',
                     'fontpad':0.0025,
                     'xborderpad':0.07,
                     'yborderpad':0.015}


m = pptmaps.PlotPhaseTensorMaps(fn_list = elst,
                                plot_freq = plot_freq ,
                                arrow_legend_dict = arrow_legend_dict,
                                ftol = 0.2,
                                xpad = 0.02,
                                plot_tipper = 'yr',
                                arrow_dict = arrow_dict,
                                ellipse_dict = ellipse_dict,
                                figsize=(8,4))

if save:
    plt.savefig(savepath,dpi=300)
