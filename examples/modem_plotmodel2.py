# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author:    u64125
@author:    Fei.zhang@ga.gov.au
"""

import os.path as op

from mtpy.imaging.plot_depth_slice import PlotDepthSlice
from mtpy.imaging.plot_phase_tensor_map import PlotPTMaps
from mtpy.imaging.plot_response import PlotResponse
from mtpy.imaging.plot_rms_map import PlotRMSMaps


# plot depth slice
# distance from slice to grid centre locations

def main(data_dir, plot_type='PTMap', di=20, periodin=0):
    # Alison wd = r'V:\Software\mtpy\development\modem_plotting\VicSynthetic07'

    # workdir = r'E:\Githubz\mtpy2\examples\data\ModEM_files'
    # workdir = r'/Softlab/Githubz/mtpy2/examples/data/ModEM_files'
    # workdir = r'/g/data/ha3/fxz547/Githubz/mtpy2/examples/data/ModEM_files'
    # wd = op.join(workdir,'VicSynthetic07')

    wd = data_dir
    plot_type = plot_type
    di = di

    # plot phase tensor map with residuals:
    # this will NOT work, an empty figure. plt.savefig(op.join(wd,'ptmaps.png'),dpi=300,ellipse_size=40)
    if plot_type == 'PTMap':
        ptmObj = PlotPTMaps(data_fn=op.join(wd, 'ModEM_Data_noise10inv.dat'),
                            resp_fn=op.join(wd, 'Modular_MPI_NLCG_019.dat'),
                            ellipse_size=30)

        outfn=op.join(wd,'ptmaps.png')
        ptmObj.plot(period=periodin, save2file=outfn)


    # plot map of RMS values
    # python examples/modem_plotmodel2.py examples/data/ModEM_files/VicSynthetic07 RMSMap
    if plot_type == 'RMSMap':
        resfile=op.join(wd, 'Modular_MPI_NLCG_019.res')
        prmsObj = PlotRMSMaps(residual_fn=resfile,xminorticks=50000, yminorticks=50000)
                   # ,depth_index=di, save_plots='y') # these are not in func args

#        prmsObj.plot_loop(fig_format="png" )    #plot all periods and save figure

    # plot responses at a station
    #FZ: how to determine this param plot_type=['VIC029']
    #Change for consistency to the save_plots='y'  instead of the implemented: plt.savefig(op.join(wd,'response.png'),dpi=300)
    if plot_type == 'Response':
        outfn=op.join(wd,'response.png')
        pltObj = PlotResponse(data_fn=op.join(wd, 'ModEM_Data_noise10inv.dat'), plot_type=['VIC099']) #, save_plots='y')
        pltObj.plot(outfn)

    # plot depth slice
    if plot_type == 'DepthSlice':
        print("plot type is", plot_type)
        modrho = op.join(wd, 'Modular_MPI_NLCG_019.rho')
        print(modrho)

        #pltObj= PlotDepthSlice(model_fn=modrho, xminorticks=100000, yminorticks=100000, depth_index=di, save_plots='y')
        pltObj = PlotDepthSlice(model_fn=modrho, save_plots='y', depth_index=20)

        pltObj.plot()

    return


#########################################################################
# plot_type=[ PTMap RMSMap Response DepthSlice ]
# How2Run:
# python examples/modem_plotmodel2.py ./examples/data/ModEM_files/VicSynthetic07 PTMap pindex
# python examples/modem_plotmodel2.py ./examples/data/ModEM_files/VicSynthetic07 RMSMap Response DepthSlice
# ---------------------------------------
if __name__ == '__main__':

    import sys

    if len(sys.argv) <= 2:
        print("USAGE example:")
        print("python %s examples/data/ModEM_files/VicSynthetic07 [PTMap|RMSMap|Response|DepthSlice]" % (sys.argv[0]))
    elif len(sys.argv) == 3:
        data_dir = sys.argv[1]
        plot_type = sys.argv[2]

        if(plot_type not in ['PTMap','RMSMap', 'Response', 'DepthSlice']):
            print("Input Parameter plot type must be in:", ['PTMap','RMSMap', 'Response', 'DepthSlice'])

        main(data_dir, plot_type=plot_type)
    else:
        data_dir = sys.argv[1]
        plot_type = sys.argv[2]
        period_index =int(sys.argv[3])

        if (plot_type not in ['PTMap', 'RMSMap', 'Response', 'DepthSlice']):
            print("Input Parameter plot type must be in:", ['PTMap', 'RMSMap', 'Response', 'DepthSlice'])

        main(data_dir, plot_type=plot_type, periodin=period_index)

