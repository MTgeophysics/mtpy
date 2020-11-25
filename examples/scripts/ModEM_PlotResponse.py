# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot data and responses from ModEM model.
To plot data only, comment out resp_fn line in PlotResponse


"""
import os

# from mtpy.imaging.plot_response import PlotResponse
from mtpy.modeling.modem import PlotResponse


#### Inputs ####
wd = r"C:\mtpywin\mtpy\examples\model_files\ModEM_2"
# wd = r'E:\Githubz\mtpy\examples\model_files\ModEM_2'
savepath = r"C:/tmp"

filestem = "Modular_MPI_NLCG_004"
datafn = "ModEM_Data.dat"

plot_z = False
respfn = filestem + ".dat"
station = "Synth02"

ro = PlotResponse(
    data_fn=os.path.join(wd, datafn),
    resp_fn=os.path.join(wd, respfn),
    plot_type=[station],
    plot_style=3,  # 1 for 4-columns; 2 for 2-columns, 3 for
    # 1-column with diagonals semi-transparent
    plot_z=False,
    res_limits=[1e-2, 1e4],
    phase_limits=[-180, 180],
    shift_yx_phase=False,
    fig_size=[3, 6],
    font_size=10,
)


ro.plot()

# ro.save_figure(os.path.join(savepath,'response.png'),
#               fig_dpi=400) # change fig_dpi to your desired resolution
