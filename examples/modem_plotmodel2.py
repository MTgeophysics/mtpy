# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: u64125
"""

import os.path as op

#refactored import mtpy.modeling.modem_new as mtmn

from mtpy.modeling.plot_phase_tensor_map import PlotPTMaps

from mtpy.modeling.plot_rms_map import PlotRMSMaps

from mtpy.modeling.plot_response import PlotResponse

from mtpy.modeling.plot_depth_slice import PlotDepthSlice


# plot depth slice
# distance from slice to grid centre locations

def main(data_dir, plot_type='PTMap', di=20):

    #Alison wd = r'V:\Software\mtpy\development\modem_plotting\VicSynthetic07'

    # workdir = r'E:\Githubz\mtpy2\examples\data\ModEM_files'
    # workdir = r'/Softlab/Githubz/mtpy2/examples/data/ModEM_files'
    # workdir = r'/g/data/ha3/fxz547/Githubz/mtpy2/examples/data/ModEM_files'
    # wd = op.join(workdir,'VicSynthetic07')

    wd=data_dir
    plot_type = plot_type
    di = di

    # plot phase tensor map with residuals
    if plot_type == 'PTMap':
        PlotPTMaps(data_fn=op.join(wd, 'ModEM_Data_noise10inv.dat'),
                        resp_fn=op.join(wd, 'Modular_MPI_NLCG_019.dat')
                        )
        # this will save an empty figureL plt.savefig(op.join(wd,'ptmaps.png'),dpi=300,ellipse_size=40)

    # plot map of RMS values
    if plot_type == 'RMSMap':
        PlotRMSMaps(residual_fn=op.join(wd, 'Modular_MPI_NLCG_019.res'),
                           xminorticks=50000, yminorticks=50000, depth_index=di, save_plots='y')

    # plot responses at a station
    if plot_type == 'Response':
        PlotResponse(data_fn=op.join(wd,'ModEM_Data_noise10inv.dat'),plot_type=['VIC029'] ,save_plots='y')
        
        #use the save_plots='y'  instead of this one: plt.savefig(op.join(wd,'response.png'),dpi=300)

    # plot depth slice
    if plot_type == 'DepthSlice':
        print("plot type is", plot_type)
        modrho=op.join(wd,'Modular_MPI_NLCG_019.rho')
        print(modrho)

        PlotDepthSlice(model_fn=modrho, xminorticks=50000,yminorticks=50000,depth_index=di,save_plots='y')


    return

#########################################################################
#plot_type=[ PTMap RMSMap Response DepthSlice ]
# How2Run:
#python examples/modem_plotmodel2.py ./examples/data/ModEM_files/VicSynthetic07 PTMap RMSMap Response DepthSlice
#---------------------------------------
if __name__ == '__main__':

    import sys

    if len(sys.argv)<2:
        print("USAGE example: python %s /g/data/ha3/fxz547/Githubz/mtpy2/examples/data/ModEM_files/VicSynthetic07"%(sys.argv[0]))
    elif len(sys.argv)==2:
        data_dir=sys.argv[1]
        main(data_dir, plot_type='PTMap')
    else:
        data_dir=sys.argv[1]
        plot_type=sys.argv[2]
        main(data_dir, plot_type=plot_type)
