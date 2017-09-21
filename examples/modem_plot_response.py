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
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

# import section
import os
# import matplotlib.pyplot as plt
# import numpy as np
# from mtpy.modeling.modem_data import Data
from mtpy.imaging.plot_response import PlotResponse

def plot_response():
    #### Default Inputs ####
    modem_data_dir = r'E:\Githubz\ModEM_plotResponse_Issue\ModEM_files'
    filestem = 'Modular_MPI_NLCG_100'
    datafn = 'ModEM_Data.dat'
    station_list = ['GB%02i'%n for n in xrange(1,40)]     # ['GB01', 'GB02',....,'GB39']
    plot_z = False

    respfn = filestem+'.dat'

    for station in station_list[0:1]:

        # plot responses at a station
        robj = PlotResponse(data_fn=os.path.join(modem_data_dir, datafn),
                          resp_fn=os.path.join(modem_data_dir, respfn),  #filestem+'.dat'),
                          plot_type=[station],
                          plot_z=plot_z,
                          #  ctmm='r',ctem='b',
                          res_limits=(.01,1000)
                          )

        # save to different image formats to compare their quality. JPG is best
        robj.plot(save2file=r'E:/tmp/test_plot_resp/plot_with_res_limits.jpg')
        robj.plot(save2file=r'E:/tmp/test_plot_resp/plot_with_res_limits.png')
        robj.plot(save2file=r'E:/tmp/test_plot_resp/plot_with_res_limits.eps')
        robj.plot(save2file=r'E:/tmp/test_plot_resp/plot_with_res_limits.pdf')
        robj.plot(save2file=r'E:/tmp/test_plot_resp/plot_with_res_limits.svg') # will create a big html file

def main():
    """
    define my main function
    :return:
    """
    print("Template main()")

    return


# =============================================
# Section for quick test of this script
# ---------------------------------------------
if __name__ == "__main__":
    # call functions
    # main()

    plot_response()
