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

# import legacy.plotptmaps as pptmaps
import mtpy.imaging.phase_tensor_maps as pptmaps
from mtpy.core.mt import MT
from tests import EDI_DATA_DIR, EDI_DATA_DIR2
from tests.imaging import ImageTestCase


class test_plotPhaseTensorMaps(ImageTestCase):
    def test_edi_files(self):
        """
        test fun
        :return:
        """
        # directory containing edis
        edipath = EDI_DATA_DIR
        # whether or not to save the figure to file
        save = True

        # full path to file to save to
        savepath = os.path.join(self._temp_dir, "phase_tensor_map.png")

        # frequency to plot
        plot_freq = 1e-2

        # gets edi file names as a list
        elst = [op.join(edipath, f) for f in os.listdir(edipath) if f.endswith(".edi")]
        mtlist = [MT(ff) for ff in elst]

        # parameters describing ellipses
        ellipse_dict = {
            "ellipse_size": 0.01,
            "ellipse_colorby": "phimax",
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
            mt_object_list=mtlist,
            plot_freq=plot_freq,
            # ftol = .5,
            # xpad = 0.02,
            # plot_tipper = 'yr',
            # ellipse_size=ellipse_dict['ellipse_size'],
            # ellipse_dict=ellipse_dict, # old line ( 22/02/2018 )
            ellipse_size=0.01,
            ellipse_colorby="phimax",
            ellipse_range=(0, 90, 1),
            ellipse_cmap="mt_bl2gr2rd",
        )
        # need to set properties and redraw
        phase_tensor_map.ellipse_size = 0.01
        phase_tensor_map.redraw_plot()

        if save:
            phase_tensor_map.save_figure(savepath, close_plot="n")
            assert os.path.exists(savepath)

    def test_edi_files2(self):
        """
        test fun
        :return:
        """

        # directory containing edis
        edipath = EDI_DATA_DIR2
        # whether or not to save the figure to file
        save = True

        # full path to file to save to
        savepath = os.path.join(self._temp_dir, "phase_tensor_map_2.png")

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
            #  ftol = .5,
            #  xpad = 0.02,
            plot_tipper="yr",
            # arrow_dict=arrow_dict, # Temporary Commented for testing
            # ellipse_dict=ellipse_dict, #Temporary Commented for testing
            #           # New incluseion ( 22/02/2018 )
            # ellipse parameters
            ellipse_size=0.1,
            ellipse_colorby="phimin",
            ellipse_range=(0, 90, 1),
            ellipse_cmap="mt_bl2gr2rd",
            # arrow parameters
            arrow_size=0.02,
            arrow_lw=0.01,
            arrow_head_width=0.002,
            arrow_head_length=0.002,
            arrow_color_real="b",
            arrow_direction=0,
            arrow_threshold=0.8,
        )

        phase_tensor_map.ellipse_size = 0.5
        phase_tensor_map.arrow_size = 10
        phase_tensor_map.redraw_plot()
        if save:
            phase_tensor_map.save_figure(savepath, close_plot="n")
            assert os.path.exists(savepath)
