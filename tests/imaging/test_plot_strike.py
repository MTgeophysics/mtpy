# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots strike
"""
import os

from mtpy.imaging.plotstrike import PlotStrike
import os.path as op

# path to edis
from tests import EDI_DATA_DIR
from tests.imaging import ImageTestCase, plt_wait, plt_close


class Test_PlotStrike(ImageTestCase):
    def test_edi_files(self):
        epath = EDI_DATA_DIR

        elst = [
            op.join(epath, edi) for edi in os.listdir(epath) if edi.endswith(".edi")
        ][::4]

        plotstrike = PlotStrike(fn_list=elst)
