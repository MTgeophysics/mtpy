#!/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    plot the responses curves from ModEM output data.
    
References: 
    https://gajira.atlassian.net/browse/ALAMP-77

CreationDate:   22/09/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     22/09/2017   FZ
    LastUpdate:     01/12/2017   FZ    debug latitude value error. Had to swap the > lat lon pair in the dat file
"""

import os

# from legacy.plot_response import PlotResponse
from mtpy.modeling.modem.plot_response import PlotResponse


def plot_response():
    #### Default Inputs ####
    modem_data_dir = r'E:\Githubz\ModEM_plotResponse_Issue\ModEM_files'
    filestem = 'Modular_MPI_NLCG_100'
    modem_data_dir = r'E:\Githubz\example_plot_response'
    filestem = 'Modular_MPI_NLCG_094.dat'
    datafn = 'ModEM_Data.dat'
    station_list = ['GB%02i' % n for n in xrange(1, 40)]  # ['GB01', 'GB02',....,'GB39']
    plot_z = False

    respfn = filestem + '.dat'

    #for station in station_list[8:10]:
    for station in ['GB08','GB09']:

        # plot responses at a station

        resp_range = None
        # resp_range = (0.01, 10000)  # This limit should be big enough, otherwise the plot curve will be out.
        if resp_range is None:
            outfile = r'./temp/plot_responses_NO_yrange.jpg'
        else:
            outfile = r'./temp/plot_responses_with_yrange.jpg'

        robj = PlotResponse(data_fn=os.path.join(modem_data_dir, datafn),
                            resp_fn=os.path.join(modem_data_dir, filestem),
                            plot_type=[station],
                            plot_style=2,
                            plot_z=plot_z,
                            #  ctmm='r',ctem='b',
                            res_limits=resp_range
                            )


        #robj.plot()

        robj.plot()

# =============================================
# Section for quick test of this script
# python examples/cmdline/modem_plot_response.py
# ---------------------------------------------
if __name__ == "__main__":
    # call functions
    # main()

    plot_response()
