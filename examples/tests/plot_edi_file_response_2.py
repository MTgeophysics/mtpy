# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots edi files (res/phase vs period) for all edis in a directory and saves out as png files
"""

import os
os.chdir(r'C:\Git\mtpy')

import mtpy.imaging.plotresponse as mtpr
import mtpy.core.edi as mtedi
import os.path as op


# path to edis
epath = r'C:\Git\mtpy\examples\data\edi_files_2'

svdir = r'C:\Git\mtpy\examples\plots\edi_plots'

elst=[op.join(epath,edi) for edi in os.listdir(epath) if (edi.endswith('.edi'))]


for efile in elst[-1:]:
    eo = mtedi.Edi(efile)
    pr = mtpr.PlotResponse(fn=efile,
                           plot_num=2,
                           plot_tipper='yri',
                           plot_pt='y')
    pr.save_plot(op.join(svdir,op.join(svdir,op.basename(efile)[:-4]+'.png')))