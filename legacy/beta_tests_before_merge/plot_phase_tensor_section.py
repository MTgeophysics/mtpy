# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period) 

fails with error:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 714, in runfile
    execfile(filename, namespace)
  File "C:\Users\u64125\AppData\Local\Continuum\Miniconda2\envs\mtpy27\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 74, in execfile
    exec(compile(scripttext, filename, 'exec'), glob, loc)
  File "U:/Software/mtpy/development/plot_phase_tensor_section.py", line 33, in <module>
    dpi=300)
  File "mtpy\imaging\plotptpseudosection.py", line 345, in __init__
    self._read_ellipse_dict()
TypeError: _read_ellipse_dict() takes exactly 2 arguments (1 given)

"""
import os
os.chdir(r'C:\Git\mtpy')
import mtpy.imaging.plotptpseudosection as ptp
import os.path as op

# path to edis
epath = r'C:\Git\mtpy\examples\data\edi_files'


elst=[op.join(epath,edi) for edi in os.listdir(epath) if edi.endswith('.edi')]

ptp.PlotPhaseTensorPseudoSection(
#mt_object_list=mtlist,
                                 fn_list = elst,
                                 tscale = 'period',
                                 ylim = (1e-1,1e3),
                                 stretch=(2,1), # determines (x,y) aspect ratio of plot
                                 station_id=(0,10), # indices for showing station names
#                                 ellipse_dict={'ellipse_size':0.5,'ellipse_colorby':'skew_seg','ellipse_range':(-12,12,3)},#,'colorby':'skew_seg','range':(-12,12,3)
                                 plot_tipper = 'yr',
                                 arrow_dict = {'size':3,'head_length':0.1,
                                               'head_width':0.1,'lw':0.5},# arrow parameters, adjust as necessary. lw = linewidth
                                 font_size=4,
                                 dpi=300)