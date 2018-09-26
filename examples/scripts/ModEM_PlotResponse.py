# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot data and responses from ModEM model.
To plot data only, comment out resp_fn line in PlotResponse


"""
import os.path as op

#from mtpy.imaging.plot_response import PlotResponse
from mtpy.modeling.modem import PlotResponse


#### Inputs ####
wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'
#wd = r'E:\Githubz\mtpy\examples\model_files\ModEM_2'
savepath = r'C:/tmp'

filestem = 'Modular_MPI_NLCG_004'
datafn = 'ModEM_Data.dat'

plot_z = False
respfn = filestem+'.dat'
station = 'Synth02'

ro = PlotResponse(data_fn=op.join(wd,datafn),
                  resp_fn=op.join(wd,respfn),
                  plot_type=[station],
                  plot_style=2,  # 1 for 4-colums; 2 for 2-columns
                  plot_z=False,
#                 fig_size=[3,2],
#                 font_size=4
                  )


ro.plot()

ro.save_figure(r'U:\Software\mtpy\example_plots'
               ,fig_dpi=400) # change fig_dpi to your desired resolution
