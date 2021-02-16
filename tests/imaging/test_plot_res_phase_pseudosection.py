# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots resistivity and phase as a coloured pseudo section (distance along profile vs period) 
"""
import os.path as op

import os

from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection
from tests import EDI_DATA_DIR, EDI_DATA_DIR2

# import matplotlib.pyplot as plt
# plt.ion() # make figure disappear automatically:
# plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc
from tests.imaging import ImageTestCase


class Test_plotResPhasePseudoSection(ImageTestCase):
    def test_edi_files(self):
        """
        test fun
        :return:
        """

        # path to edis
        epath = EDI_DATA_DIR

        save_path = os.path.join(self._temp_dir, "resphase.png")

        elst = [
            op.join(epath, edi) for edi in os.listdir(epath) if edi.endswith(".edi")
        ][::4]

        print(elst)
        resphase = PlotResPhasePseudoSection(fn_list=elst)

        resphase.save_plot(save_path, close_plot="n")

        assert os.path.exists(save_path)

    def test_edi_files2(self):
        # path to edis
        epath = EDI_DATA_DIR2

        save_path = os.path.join(self._temp_dir, "resphase_2.png")

        elst = [
            op.join(epath, edi) for edi in os.listdir(epath) if edi.endswith(".edi")
        ][::4]

        resphase = PlotResPhasePseudoSection(fn_list=elst)

        resphase.save_plot(save_path, close_plot="n")

        assert os.path.exists(save_path)
