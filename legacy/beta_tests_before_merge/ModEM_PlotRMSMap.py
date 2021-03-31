# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot RMS at each station as a map

"""
import os

os.chdir(r"C:\Git\mtpy")
from mtpy.modeling.modem import PlotRMSMaps

# directory where files are located
wd = r"C:\Git\mtpy\examples\model_files\ModEM"

# directory to save to
save_path = r"C:\Git\mtpy\examples\plots\ModEM"

# file stem for inversion result
filestem = "Modular_MPI_NLCG_004"

# period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
period_index = 0

# plot map
rmsmap = PlotRMSMaps(
    residual_fn=os.path.join(wd, filestem + ".res"),
    period_index=period_index,
    xminorticks=50000,
    yminorticks=50000,
    save_plots="y",
)
rmsmap.save_figure(save_path)
