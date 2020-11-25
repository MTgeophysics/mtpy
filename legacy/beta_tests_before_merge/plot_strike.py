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
  File "C:/Git/mtpy/examples/tests/plot_strike.py", line 21, in <module>
    plotstrike = PlotStrike(fn_list=elst)
  File "mtpy\imaging\plotstrike.py", line 240, in __init__
    self.plot()
  File "mtpy\imaging\plotstrike.py", line 307, in plot
    zinv = mt.get_Zinvariants()
AttributeError: 'MTplot' object has no attribute 'get_Zinvariants'

"""
import os

os.chdir(r"C:\Git\mtpy")
from mtpy.imaging.plotstrike import PlotStrike
import os.path as op
import matplotlib.pyplot as plt

# path to edis
epath = r"C:\Git\mtpy\examples\data\edi_files"


elst = [op.join(epath, edi) for edi in os.listdir(epath) if edi.endswith(".edi")][::4]

plotstrike = PlotStrike(fn_list=elst)
