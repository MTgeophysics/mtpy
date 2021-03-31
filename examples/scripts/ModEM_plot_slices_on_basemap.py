# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:53:46 2019

@author: u64125
"""

from mtpy.modeling.modem import PlotSlices
from tests import MODEM_DIR, TEST_TEMP_DIR

model_fn = MODEM_DIR.joinpath("Modular_MPI_NLCG_004.rho")
data_fn = MODEM_DIR.joinpath("ModEM_Data.dat")


ps = PlotSlices(
    model_fn,
    data_fn=data_fn,
    model_epsg=28353,  # epsg used in projection of model grid to ensure correct projection back to lat/long grid
    cmap="jet_r",  # color map
    climits=(0, 4),  # log10(color_limits) for resistivity
    plot_stations=True,  # True/False, whether or not to plot stations
    draw_colorbar=True,  # whether or not to show a colorbar
    plot_yn="n",
)

# loop through depths
for depth in [2e3, 10e3, 20e3]:
    ps.basemap_plot(
        depth,  # depth to plot (code finds nearest slice)
        tick_interval=None,  # tick interval in degrees, if not provided or None it is calculated from data
        buffer=None,  # buffer around stations in degrees, if not provided or None it is calculated from data
        mesh_rotation_angle=0,  # option to specify the mesh rotation angle, if rotated grid was used
        save=True,
        save_path=TEST_TEMP_DIR,  # savepath to save figure to. If not provided or None/False, figure is not saved
    )
