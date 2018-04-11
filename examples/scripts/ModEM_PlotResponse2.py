# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot data and responses from ModEM model.
To plot data only, comment out resp_fn line in PlotResponse


"""
import os.path as op
import os
#os.chdir(r'C:/mtpywin/mtpy')
#from mtpy.imaging.plot_response import PlotResponse
from mtpy.modeling.modem import PlotResponse


#### Inputs ####
wd = r'U:\RegionalSurvey\MT046_GeorginaArunta\Modelling\ModEM\GBinv52tip2'
wd = r'E:\Githubz\mtpy\zprivate\modem_plot_response_issues\GBinv52tip2'
savepath = wd# r'U:\Software\mtpy\example_plots'

filestem = 'Modular_MPI_NLCG_108'
datafn = 'ModEM_Data.dat'

plot_z = False
respfn = filestem+'.dat'

#station = 'GB08' # no tipper data, was not working for plot() with style=1 which is intended to produce 4-column figure
station = 'GB09'  # has tipper data

ro = PlotResponse(data_fn=op.join(wd,datafn),
                  resp_fn=op.join(wd,respfn),
                  plot_type=[station],
                  plot_style=1,  # 1 for 4-colums; 2 for 2-columns
                  plot_z=False,
                 # fig_size=[3,2],
                 # font_size=4
                  )


# print("calling  plot() ....")
ro.plot()

#ro.save_figure(op.join(wd,))
