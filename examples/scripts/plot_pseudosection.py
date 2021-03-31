# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period) 
"""
import os.path as op
import os

os.chdir(r"C:\mtpywin\mtpy")

import mtpy.imaging.plotpseudosection as pps


epath = r"C:\mtpywin\mtpy\examples\data\edi_files"

elst = [
    op.join(epath, edi) for edi in os.listdir(epath) if (edi.endswith(".edi"))
]  # and edi.startswith('GB')

pps.PlotResPhasePseudoSection(
    fn_list=elst,
    linedir="ns",
    plot_xx="n",
    plot_xy="y",
    plot_yx="y",
    plot_yy="n",
    res_limits=[0, 3],
    phase_limits=[0, 90],
    shift_yx_phase=True,
    plot_style="pcolormesh",
)
