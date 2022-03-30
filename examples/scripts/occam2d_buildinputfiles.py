# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:51:02 2015

@author: Alison Kirkby

sets up input files for running 2d occam inversions using the occam2d_rewrite module

"""

import mtpy.modeling.occam2d as occam2d
import os
import os.path as op
import numpy as np


# path where edi files are located
edipath = r"C:\mtpywin\mtpy\examples\data\edi_files"

# path to save to
savepath = r"C:/tmp"

if not op.exists(savepath):
    os.mkdir(savepath)


# list of stations
slst = [edi[0:-4] for edi in os.listdir(edipath) if edi.find(".edi") > 0]


# create an occam data object
ocd = occam2d.Data(
    edi_path=edipath,
    station_list=slst,
    interpolate_freq=True,
    freq=np.logspace(-3, 3, 37),
    model_mode="log_all",
)
ocd.save_path = savepath
ocd.freq_num = 50  # number of frequencies to invert for

#### make data file
# geoelectric strike for rotation
# if not specified will calculate from the data
ocd.geoelectric_strike = 1

# error floors
ocd.res_te_err = 10
ocd.res_tm_err = 10
ocd.phase_te_err = 5
ocd.phase_tm_err = 5
# ocd.model_mode= 4
ocd.write_data_file()


# make model and mesh files
ocr = occam2d.Regularization(ocd.station_locations)
# number of layers
ocr.n_layers = 60

ocr.cell_width = 500  # cell width to aim for, note
# this is the mesh size (2 mesh
# blocks per model block)
ocr.x_pad_multiplier = 1.9  # controls size of padding
ocr.trigger = 0.25  # controls aspect ratio of blocks
# ocr.z_bottom = 200000

# z1 layer and target depth in metres
ocr.z1_layer = 50
ocr.z_target_depth = 80000
ocr.save_path = ocd.save_path
ocr.build_mesh()
ocr.build_regularization()
ocr.write_mesh_file()
ocr.write_regularization_file()
ocr.plot_mesh()

# make startup file
ocs = occam2d.Startup()
ocs.iterations_to_run = 40
ocs.data_fn = op.join(ocd.save_path, "OccamDataFile.dat")
ocs.resistivity_start = 2.0
ocr.get_num_free_params()
ocs.param_count = ocr.num_free_param
ocs.save_path = ocd.save_path
ocs.model_fn = ocr.reg_fn
ocs.write_startup_file()
