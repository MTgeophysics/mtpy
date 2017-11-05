# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 07:29:58 2013

@author: Alison Kirkby

plots phase tensor ellipses as a map for a given frequency

bug in setting ellipse properties: 
The ellipse properties are not being set via the arguments - need to create a
phase_tensor_map object then set the properties then run redraw_plot.

"""
import os
import os.path as op

import mtpy.imaging.plotptmaps as pptmaps
from mtpy.core.mt import MT
from tests.beta import *
# import matplotlib.pyplot as plt
# plt.ion()  # make figure disappear automatically:
# plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc

def test_func():
    """
    test fun
    :return:
    """
    # directory containing edis
    edipath =EDI_DATA_DIR
    # whether or not to save the figure to file
    save = True

    # full path to file to save to
    savepath = os.path.join(TEMP_OUT_DIR,'phase_tensor_map.png')

    # frequency to plot
    plot_freq = 1e-2

    # gets edi file names as a list
    elst = [op.join(edipath, f) for f in os.listdir(edipath) if f.endswith('.edi')]
    mtlist = [MT(ff) for ff in elst]

    # parameters describing ellipses
    ellipse_dict = {'ellipse_size': .01, 'ellipse_colorby': 'phimax', 'ellipse_range': (0, 90, 1),
                    'cmap': 'mt_bl2gr2rd'}

    # parameters describing the induction vector arrows
    arrow_dict = {'arrow_size': 0.02,
                  'arrow_lw': 0.01,
                  'arrow_head_width': 0.002,
                  'arrow_head_length': 0.002,
                  'arrow_color_real': 'b',
                  'direction': 0,
                  'threshold': 0.8}

    phase_tensor_map = pptmaps.PlotPhaseTensorMaps(
        mt_object_list=mtlist,
        plot_freq=plot_freq,
        # ftol = .5,
        # xpad = 0.02,
        # plot_tipper = 'yr',
        # arrow_dict = arrow_dict,
        ellipse_size=ellipse_dict['ellipse_size'],
        ellipse_dict=ellipse_dict,
    )
    # need to set properties and redraw
    phase_tensor_map.ellipse_size = 0.01
    phase_tensor_map.redraw_plot()

    if save:
        phase_tensor_map.save_figure(savepath)
        assert (os.path.exists(savepath))

#######################################
if __name__ == "__main__":
    test_func()
