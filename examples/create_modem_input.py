# -*- coding: utf-8 -*-
"""
Create modem input files:
This script includes topography in the model. To not include topography,
set number of air layers to zero (recommended) or comment out add_topography
line. Note: setting number of air layers to zero will add bathymetry but not
topography.

USAGE examples:
python examples/create_modem_input.py tests/data/edifiles/ examples/etopo1.asc /e/tmp/modem_test
python examples/create_modem_input.py /e/Data/MT_Datasets/WenPingJiang_EDI /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc
       /e/tmp/WenPingTest
python examples/create_modem_input.py /e/Data/MT_Datasets/concurry_EDI_files/ /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc
    /e/tmp/Concurry

Developed by
    Alison.Kirkby@ga.gov.au
    Fei.Zhang@ga.gov.au

Create Date: 2017-02-01
"""
from __future__ import print_function

import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from mtpy.core.edi_collection import EdiCollection
from mtpy.modeling.modem_covariance import Covariance
from mtpy.modeling.modem_data import Data
from mtpy.modeling.modem_model import Model

# YG: patch that changes the matplotlib behaviour
plt.ion()  # enable interactive
# plt.ioff()  # disable interactive, which will also disable this patch


def show_patcher(show_func):
    """
    patch the plt.show() if interactive is enabled to display and then close the plot after 1 second
    so plt.show() will not block the script and the figure is still visible to the user
    :param show_func:
    :return:
    """
    def new_show_func(*args, **kwargs):
        stuff = show_func(*args, **kwargs)
        # wait 1 second for the image to show on screen
        figManager = plt.gcf()
        if figManager is not None:
            canvas = figManager.canvas
            # if canvas.figure.stale:
            #     canvas.draw()
            # show(block=False)
            canvas.start_event_loop(1)  # wait time = 1
        plt.close()
        return stuff
    return new_show_func if plt.isinteractive() else show_func


# plt.show = show_patcher(plt.show)
# end of patch


def select_periods(edifiles_list):
    """
    FZ: Use edi_collection to analyse the whole set of EDI files
    :param edifiles:
    :return:
    """
    import matplotlib.pyplot as plt

    edis_obj = EdiCollection(edifiles_list)

    uniq_period_list = edis_obj.all_unique_periods  # filtered list of periods ?
    print("Unique periods", len(uniq_period_list))

    plt.hist(edis_obj.mt_periods, bins=uniq_period_list)
    # plt.hist(edis_obj.mt_periods, bins=1000)
    plt.title("Histogram with uniq_periods bins")
    plt.xlabel("Periods")
    plt.ylabel("Occurance in number of MT stations")
    plt.show()

    # 1 ASK user to input a Pmin and Pmax

    # 2 percetage stats
    # select commonly occured frequencies from all stations.
    # This could miss some slightly varied frequencies in the middle range.
    select_period_list = np.array(edis_obj.get_periods_by_stats(percentage=10.0))
    print("Selected periods ", len(select_period_list))

    return select_period_list


if __name__ == '__main__':

    if len(sys.argv) < 4:
        print("USAGE: %s  path2edifiles path2topo.asc path2outdir" %
              sys.argv[0])
        sys.exit(1)
    else:
        edipath = sys.argv[1]  # edi files to be inversioned
        topofile = sys.argv[2]  # topography file, if using
        outputdir = sys.argv[3]  # path to save to

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    # epsg to project to. Google epsg 'your projection'
    epsg_code = 28354
    epsg_code = 3112

    edi_list = glob.glob(edipath + '/*.edi')

    if edi_list is None or (edi_list) < 1:
        print("Error: No edi files found in the dir %s" % edipath)
        sys.exit(2)

    # period list (can take periods from one of the edi files, or just specify
    # periods directly using the logspace function (commented out))

    # eo = mtedi.Edi(edi_list[0])  # this may miss some periods?
    # period_list = 1. / eo.Z.freq # period_list = np.logspace(-3,3)

    period_list = select_periods(edi_list)

    datob = Data(edi_list=edi_list,
                 inv_mode='1',
                 period_list=period_list,
                 epsg=epsg_code,
                 error_type='floor',
                 error_floor=10)
    # period_buffer=0.000001)

    datob.write_data_file(save_path=outputdir)

    # create mesh grid model object
    model = Model(Data=datob,
                  epsg=epsg_code,  # epsg
                  # cell_size_east=500, cell_size_north=500,  # concurry
                  cell_size_east=10000, cell_size_north=10000, #GA_VIC
                  # cell_size_east=1000, cell_size_north=1000, # Concurry
                  pad_north=8,  # number of padding cells in each of the north and south directions
                  pad_east=8,  # number of east and west padding cells
                  pad_z=8,  # number of vertical padding cells
                  pad_stretch_v=1.5,  # factor to increase by in padding cells (vertical)
                  pad_stretch_h=1.5,  # factor to increase by in padding cells (horizontal)
                  n_airlayers=0,  # number of air layers 0, 10, 20, depend on topo elev height
                  res_model=100,  # halfspace resistivity value for initial reference model
                  n_layers=50,  # total number of z layers, including air and pad_z
                  z1_layer=50,  # first layer thickness metres, depend
                  z_target_depth=500000)

    model.make_mesh()  # the data file will be re-write in this method. No topo elev file used yet

    model.plot_mesh()
    model.plot_mesh_xy()
    model.plot_mesh_xz()

    # write a model file and initialise a resistivity model
    model.write_model_file(save_path=outputdir)


#=========== now add topo data, with or without air layers?
    # 1) the data file will be changed in 3 columns sxi, syi and szi meters
    # 2) The covariance file will be written.
    # 3) the model file not changed?? No air layers can be seen in the .ws file.

    # add topography, define an initial resistivity model, modify and re-write the data file, define covariance mask
    # dat file will be changed and rewritten,
    # grid centre is used as the new origin of coordinate system, topo data used in the elev column.

    # model.add_topography(topofile, interp_method='nearest')  # dat file will be written again as elevation updated

    model.add_topography_2mesh(topofile, interp_method='nearest')  # dat file will be written again as elevation updated

    model.plot_topograph()  # plot the MT stations on topography elevation data

    print("*** Re-writing model file after topo data and air layers are added - will include air sea-water resistivity")
    model.write_model_file(save_path=model.save_path)
    # model.write_model_file(save_path='temp/')


    # make covariance (mask) file
    cov = Covariance(mask_arr=model.covariance_mask,
                     save_path=outputdir,
                     smoothing_east=0.3,
                     smoothing_north=0.3,
                     smoothing_z=0.3)

    cov.write_covariance_file(model_fn=model.model_fn)
