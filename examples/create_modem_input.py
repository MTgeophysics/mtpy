# -*- coding: utf-8 -*-
"""
Crreate modem input files:
This script includes topography in the model. To not include topography,
set number of air layers to zero (recommended) or comment out add_topography
line. Note: setting number of air layers to zero will add bathymetry but not
topography.

USAGE examples:
python examples/create_modem_input.py tests/data/edifiles/ examples/etopo1.asc /e/tmp/modem_test
python examples/create_modem_input.py /e/Data/MT_Datasets/WenPingJiang_EDI /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc
       /e/tmp/WenPingTest

Developed by
    Alison.Kirkby@ga.gov.au
    Fei.Zhang@ga.gov.au

Create Date: 2017-02-01
"""
from __future__ import print_function
import os
import sys
import glob
import numpy as np
import mtpy.modeling.modem as mtmn
from mtpy.core.edi_collection import EdiCollection

if __name__ == '__main__':

    if len(sys.argv) < 4:
        print ("USAGE: %s  path2edifiles path2topo.asc path2outdir" %
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
        print ("Error: No edi files found in the dir %s" % edipath)
        sys.exit(2)

    # period list (can take periods from one of the edi files, or just specify
    # periods directly using the logspace function (commented out))

    # eo = mtedi.Edi(edi_list[0])  # this may miss some periods?
    # period_list = 1. / eo.Z.freq # period_list = np.logspace(-3,3)

    # FZ: Use edi_collection to analyse the whole set of EDI files
    edis_obj = EdiCollection(edi_list)

    period_list = edis_obj.all_unique_periods  # filtered list of periods ?
    print(type(period_list))

    # select commonly occured frequencies from all stations.
    # This could miss some slightly varied frquencies in the middle range.
    period_list = np.array(edis_obj.get_periods_by_stats(percentage=10.0))
    print(type(period_list))

    datob = mtmn.Data(edi_list=edi_list,
                      inv_mode='1',
                      period_list=period_list,
                      epsg=epsg_code,
                      error_type='floor',
                      error_floor=10)
    # period_buffer=0.000001)

    datob.write_data_file(save_path=outputdir)

    # create model file
    model = mtmn.Model(Data=datob,
                       cell_size_east=2000,
                       cell_size_north=2000,
                       pad_north=7,  # number of padding cells in each of the north and south directions
                       pad_east=7,  # number of east and west padding cells
                       pad_z=6,  # number of vertical padding cells
                       pad_stretch_v=3,
                       # factor to increase by in padding cells (vertical)
                       pad_stretch_h=3,
                       # factor to increase by in padding cells (vertical)
                       n_airlayers=10,  # number of air layers
                       res_model=100,  # halfspace resistivity value for reference model
                       n_layers=80,  # total number of z layers, including air
                       z1_layer=100,  # first layer thickness
                       epsg=epsg_code,  # epsg
                       z_target_depth=120000)

    model.make_mesh()

    model.plot_mesh()

    # write a model file to initialise a resistivity model
    model.write_model_file(save_path=outputdir)

    # add topography to res model
    model.add_topography(topofile, interp_method='nearest')

    # make covariance file
    cov = mtmn.Covariance(mask_arr=model.covariance_mask,
                          save_path=outputdir,
                          smoothing_east=0.3,
                          smoothing_north=0.3,
                          smoothing_z=0.3)

    cov.write_covariance_file(model_fn=model.model_fn)
