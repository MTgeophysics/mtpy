# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots strike

fails with error:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 714, in runfile
    execfile(filename, namespace)
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 74, in execfile
    exec(compile(scripttext, filename, 'exec'), glob, loc)
  File "C:/Git/mtpy/examples/tests/fails_plot_strike.py", line 21, in <module>
    plotstrike = PlotStrike(fn_list=elst)
  File "mtpy\imaging\plotstrike.py", line 240, in __init__
    self.plot()
  File "mtpy\imaging\plotstrike.py", line 307, in plot
    zinv = mt.get_Zinvariants()
AttributeError: 'MTplot' object has no attribute 'get_Zinvariants'

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
