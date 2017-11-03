# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots resistivity and phase as a coloured pseudo section (distance along profile vs period) 
"""
import os
os.chdir(r'C:\Git\mtpy')
from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection
import os.path as op
import matplotlib.pyplot as plt

# path to edis
epath = r'C:\Git\mtpy\examples\data\edi_files_2'

save_path = r'C:\Git\mtpy\examples\plots\edi_plots\resphase_2.png'

elst=[op.join(epath,edi) for edi in os.listdir(epath) if edi.endswith('.edi')][::4]


resphase = PlotResPhasePseudoSection(fn_list=elst)

resphase.save_plot(save_path)