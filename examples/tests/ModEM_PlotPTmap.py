# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot data and responses from ModEM model.

Fails with error:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 714, in runfile
    execfile(filename, namespace)
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 74, in execfile
    exec(compile(scripttext, filename, 'exec'), glob, loc)
  File "C:/Git/mtpy/examples/tests/ModEM_PlotPTmap.py", line 40, in <module>
    ellipse_size=20
  File "mtpy\modeling\modem.py", line 5774, in __init__
    self._read_ellipse_dict()
TypeError: _read_ellipse_dict() takes exactly 2 arguments (1 given)

"""

import os
import os.path as op

os.chdir(r'C:\Git\mtpy')
from mtpy.modeling.ModEM import PlotPTMaps

wd = r'C:\Git\mtpy\examples\model_files\ModEM'

filestem = 'Modular_MPI_NLCG_004'
datafn = 'ModEM_Data.dat'


PlotPTMaps(data_fn = op.join(wd,datafn),
                resp_fn = op.join(wd,filestem + '.dat'),
                ellipse_size=20
                )
