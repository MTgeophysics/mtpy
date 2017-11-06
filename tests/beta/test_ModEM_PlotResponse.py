# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot responses from ModEM model.

fails with error:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 714, in runfile
    execfile(filename, namespace)
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 74, in execfile
    exec(compile(scripttext, filename, 'exec'), glob, loc)
  File "C:/Git/mtpy/examples/tests/ModEM_PlotResponse.py", line 49, in <module>
    plot_z=plot_z,
  File "mtpy\modeling\modem.py", line 4929, in __init__
    self.plot()
  File "mtpy\modeling\modem.py", line 5015, in plot
    z_obj._compute_res_phase()
AttributeError: 'Z' object has no attribute '_compute_res_phase'

"""
import os
import os.path as op
from unittest import TestCase

from mtpy.modeling.modem import PlotResponse

from tests.beta import SAMPLE_DIR
from tests.imaging import _plt_wait, _plt_close


class Test_ModEM_PlotResponse(TestCase):
    def tearDown(self):
        _plt_wait(5)
        _plt_close()

    def test_modular_MPI_NLCG_004(self):
        wd = op.normpath(op.join(SAMPLE_DIR, 'ModEM'))
        filestem = 'Modular_MPI_NLCG_004'
        datafn = 'ModEM_Data.dat'
        station = 'pb23'
        plot_z = False

        ro = PlotResponse(data_fn=op.join(wd, datafn),
                          resp_fn=op.join(wd, filestem + '.dat'),
                          plot_type=[station],
                          plot_z=plot_z)

        ro.plot()
