# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots edi files (res/phase vs period) for all edis in a directory and saves out as png files
"""

import os
from tests.beta import *
import mtpy.imaging.plotresponse as mtpr
import matplotlib.pyplot as plt
plt.ion() # make figure disappear automatically:
#plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc


def test_func():
    # path to edis
    epath = EDI_DATA_DIR

    svdir = TEMP_OUT_DIR

    elst=[os.path.join(epath,edi) for edi in os.listdir(epath) if (edi.endswith('.edi'))]


    for efile in elst[:3]:
        # eo = mtedi.Edi(efile)
        pr = mtpr.PlotResponse(fn=efile,
                               plot_num=2,
                               plot_tipper='yri',
                               plot_pt='y')

        figfile = os.path.join(svdir, os.path.basename(efile)[:-4]+'.png')
        pr.save_plot( figfile )

        assert (os.path.exists(figfile))


if __name__ == "__main__":
    test_func()

