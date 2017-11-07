# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots strike
"""
import os
from unittest import TestCase

from mtpy.imaging.plotstrike import PlotStrike
import os.path as op

# path to edis
from tests.beta import EDI_DATA_DIR
from tests.imaging import _plt_wait, _plt_close


class Test_PlotStrike(TestCase):
    def tearDown(self):
        _plt_wait(5)
        _plt_close()

    def test_edi_files(self):
        epath = EDI_DATA_DIR

        elst = [op.join(epath, edi) for edi in os.listdir(epath) if edi.endswith('.edi')][::4]

        plotstrike = PlotStrike(fn_list=elst)
