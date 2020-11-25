# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 07:29:58 2013

@author: Alison Kirkby

plots phase tensor ellipses as a map for a given frequency

bug in setting ellipse properties: 
The ellipse properties are not being set via the arguments - need to create a
phase_tensor_map object then set the properties then run redraw_plot.

bug/quirk in save_figure:
opens up a new figure each time it runs


"""
import os

os.chdir(r"C:\Git\mtpy")

import os.path as op
import legacy.plotptmaps as pptmaps
from mtpy.core.mt import MT

# directory containing edis
edipath = r"C:\Git\mtpy\examples\data\edi_files_2"
# whether or not to save the figure to file
save = True

# full path to file to save to
savepath = r"C:\Git\mtpy\examples\plots\edi_plots\phase_tensor_map_2.png"

# frequency to plot
plot_freq = 1.318400e-01

# gets edi file names as a list
elst = [op.join(edipath, f) for f in os.listdir(edipath) if f.endswith(".edi")]
mtlist = [MT(ff) for ff in elst]

# parameters describing ellipses
ellipse_dict = {
    "ellipse_size": 0.1,
    "ellipse_colorby": "phimin",
    "ellipse_range": (0, 90, 1),
    "cmap": "mt_bl2gr2rd",
}

# parameters describing the induction vector arrows
arrow_dict = {
    "arrow_size": 0.02,
    "arrow_lw": 0.01,
    "arrow_head_width": 0.002,
    "arrow_head_length": 0.002,
    "arrow_color_real": "b",
    "direction": 0,
    "threshold": 0.8,
}


phase_tensor_map = pptmaps.PlotPhaseTensorMaps(
    # fn_list = elst,
    mt_object_list=mtlist,
    plot_freq=plot_freq,
    #                                ftol = .5,
    #                                xpad = 0.02,
    plot_tipper="yr",
    arrow_dict=arrow_dict,
    ellipse_dict=ellipse_dict,
)

phase_tensor_map.ellipse_size = 0.5
phase_tensor_map.arrow_size = 10
phase_tensor_map.redraw_plot()
if save:
    phase_tensor_map.save_figure(savepath)
