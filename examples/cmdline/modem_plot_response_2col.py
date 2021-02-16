#! /usr/bin/env python
"""
Description:
    Alison said I found another bug in plot_response. (plot_2col function).
    Please see the attached code (GB08-1.zip) and data (data.zip).
    In the plot with tipper data, the resistivity xx and yy axes labels are incorrect.
    They should read 10-3, 10-2 etc. Both resistivity axes labels (i.e. xx yy and xy yx) resist be in scientific notation e.g. 100, 101 etc like the plot with no tipper data. Could you please fix?


    
References: 
    https://gajira.atlassian.net/browse/AUSLAMP-175

CreationDate:   24/04/2018

Developer:      fei.zhang@ga.gov.au
"""


import os
from mtpy.modeling.modem import PlotResponse


if __name__ == "__main__":

    #### Inputs ####
    wd = r"E:\Githubz\Alison_Bugs\data"
    savepath = r"E:\Githubz\Alison_Bugs\output"

    filestem = "Modular_MPI_NLCG_108"
    datafn = "ModEM_Data.dat"

    Resist_Only = False
    # True to plot impedance,
    # False for plotting resistivity and phase

    respfn = filestem + ".dat"
    station = ["GB08", "GB09"]

    ro = PlotResponse(
        data_fn=os.path.join(wd, datafn),
        resp_fn=os.path.join(wd, respfn),
        plot_type=station,
        plot_style=1,
        plot_z=Resist_Only,
        save_plots=True,
        ctem="b",
        ctmm="r",
        #                  mtem=
        #                  plot_yn = False,
        #                  fig_size=[3,2],
        #                  font_size=4
    )
    ro.plot()
