# -*- coding: utf-8 -*-
"""
Plot PTMap|RMSMap|Response|DepthSlice from Inversion Model Results

Created on Tue Oct 04 13:13:29 2016

@author:    Alison.Kirkby@ga.gov.au
@author:    Fei.zhang@ga.gov.au

Test Run:
python examples/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07
"""
import os
import sys

from mtpy.imaging.modem_phase_tensor_maps import PlotPTMaps
from mtpy.imaging.plot_depth_slice import PlotDepthSlice
from mtpy.imaging.plot_response import PlotResponse
from mtpy.imaging.plot_rms_map import PlotRMSMaps

# original test case:
# datfn='ModEM_Data_noise10inv.dat'  # what is this noiseinv.dat?
# NLCG_datfn='Modular_MPI_NLCG_019.dat'
# resfn='Modular_MPI_NLCG_019.res'
# rhofn='Modular_MPI_NLCG_019.rho'

# FZ: below works fine
# datfn='Isa_run3_NLCG_049.dat' #'ModEM_Data_noise10inv.dat'
# NLCG_datfn='Isa_run3_NLCG_049.dat'
# resfn='Isa_run3_NLCG_049.res'
# rhofn='Isa_run3_NLCG_049.rho'

# rename/copy the final MODEM results to these file names:
datfn='NLCG.dat'  # 'ModEM_Data_noise10inv.dat'
NLCG_datfn='NLCG.dat'
resfn='NLCG.res'
rhofn='NLCG.rho'


def plot_model(data_dir, plot_type='PTMap', depth_index=20, periodin=0):
    """
    plot model of the plot_type
    :param data_dir: directory where modem's NLCG.dat .rho .res files are located
    :param plot_type: one of these 4: PTMap|RMSMap|Response|DepthSlice
    :param di:
    :param periodin:
    :return:
    """

    wd=data_dir
    plot_type=plot_type
    depth_index = depth_index  # depth index

    # plot phase tensor map with residuals:
    # this will NOT work, an empty figure.
    # plt.savefig(op.join(wd,'ptmaps.png'),dpi=300,ellipse_size=40)
    if plot_type == 'PTMap':
        ptmObj=PlotPTMaps(data_fn=os.path.join(wd, datfn),
                          resp_fn=os.path.join(wd, NLCG_datfn),
                          ellipse_size=30)

        outfn=os.path.join(wd, 'ptmaps.png')
        ptmObj.plot(period=periodin, save2file=outfn)

    # plot map of RMS values
    # python examples/modem_plotmodel2.py
    # examples/data/ModEM_files/VicSynthetic07 RMSMap
    if plot_type == 'RMSMap':
        resfile=os.path.join(wd, resfn)
        prmsObj=PlotRMSMaps(
            residual_fn=resfile,
            xminorticks=50000,
            yminorticks=50000)
        # ,depth_index=di, save_plots='y') # these are not in func args

    # prmsObj.plot_loop(fig_format="png" )    #plot all periods and save
    # figure

    # plot responses at a station
    if plot_type == 'Response':
        outfn = os.path.join(wd, 'response.png')
        pltObj=PlotResponse(data_fn=os.path.join(wd, datfn),plot_type=['16-L03S01','VIC001'])
        #FZ: need to refactor plot_type= list of  station names

        pltObj.plot(outfn)

    # plot depth slice
    if plot_type == 'DepthSlice':
        print("plot type is", plot_type)
        modrho=os.path.join(wd, rhofn)
        print(modrho)

        # pltObj= PlotDepthSlice(model_fn=modrho, xminorticks=100000, yminorticks=100000, depth_index=di, save_plots='y')
        pltObj=PlotDepthSlice(
            model_fn=modrho,
            save_plots='y',
            depth_index=depth_index)

        pltObj.plot(ind=depth_index)

    return


#########################################################################
# plot_type=[ PTMap RMSMap Response DepthSlice ]
# How2Run:
# python examples/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07 PTMap pindex
# python examples/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07 RMSMap Response DepthSlice
# ---------------------------------------
if __name__ == '__main__':

    if len(sys.argv) <= 2:
        print("USAGE example:")
        print(
            "python %s examples/data/ModEM_files/VicSynthetic07 [PTMap|RMSMap|Response|DepthSlice]" %
            (sys.argv[0]))
        for plot_type in ['PTMap', 'RMSMap', 'Response', 'DepthSlice']:
            plot_model(sys.argv[1], plot_type=plot_type)

    elif len(sys.argv) == 3:
        data_dir=sys.argv[1]
        plot_type=sys.argv[2]

        if (plot_type not in ['PTMap', 'RMSMap', 'Response', 'DepthSlice']):
            print("Input Parameter plot type must be in:", [
                'PTMap', 'RMSMap', 'Response', 'DepthSlice'])

        plot_model(data_dir, plot_type=plot_type)
    else:
        data_dir=sys.argv[1]
        plot_type=sys.argv[2]
        period_index=int(sys.argv[3])

        if (plot_type not in ['PTMap', 'RMSMap', 'Response', 'DepthSlice']):
            print("Input Parameter plot type must be in:", [
                'PTMap', 'RMSMap', 'Response', 'DepthSlice'])

        plot_model(data_dir, plot_type=plot_type, periodin=period_index)
