# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period)
"""
import os
import os.path as op

import mtpy.imaging.phase_tensor_pseudosection as ptp

# path to edis
from tests import EDI_DATA_DIR
from tests.imaging import ImageTestCase


class test_PlotPtPseudoSection(ImageTestCase):
    def test_edifiles(self):
        epath = EDI_DATA_DIR

        elst = [
            op.join(epath, edi) for edi in os.listdir(epath) if edi.endswith(".edi")
        ]

        plt_obj = ptp.PlotPhaseTensorPseudoSection(
            # mt_object_list=mtlist,
            fn_list=elst,
            tscale="period",
            ylim=(1e-1, 1e3),
            stretch=(2, 1),  # determines (x,y) aspect ratio of plot
            station_id=(0, 10),  # indices for showing station names
            #   ellipse_dict={'ellipse_size':0.5,'ellipse_colorby':'skew_seg','ellipse_range':(-12,12,3)},#,'colorby':'skew_seg','range':(-12,12,3)
            plot_tipper="yr",
            arrow_dict={
                "size": 3,
                "head_length": 0.1,
                "head_width": 0.1,
                "lw": 0.5,
            },  # arrow parameters, adjust as necessary. lw = linewidth
            font_size=4,
            dpi=300,
        )
        plt_obj.plot()
