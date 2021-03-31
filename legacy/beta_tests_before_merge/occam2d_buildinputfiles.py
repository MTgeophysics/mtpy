# -*- coding: utf-8 -*-


import os

os.chdir(r"C:/Git/mtpy")

import mtpy.modeling.occam2d as occam2d
import os.path as op
import numpy as np


# path where edi files are located
edipath = r"C:\Git\mtpy\examples\data\edi_files"

# path to save to
savepath = r"C:\Git\mtpy\examples\model_files\Occam2d"

if not op.exists(savepath):
    os.mkdir(savepath)


# list of stations
slst = [edi[0:-4] for edi in os.listdir(edipath) if edi.find(".edi") > 0]


# create an occam data object
ocd = occam2d.Data(
    edi_path=edipath,
    station_list=slst,
    #                  interpolate_freq=True,
    #                  freq=np.logspace(-3,1,30)
)


ocd.save_path = savepath

# choose frequency range to invert

# ocd.freq_num = 50
ocd.freq_min = 1
ocd.freq_max = 10000
# ocd.freq_num = 50 # number of frequencies to invert for
###########make data file

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
# cell width to aim for, note this is the mesh size (2 mesh blocks per model block)
ocr.cell_width = 200
# controls number and size of padding
ocr.num_x_pad_cells = 9
ocr.x_pad_multiplier = 1.9
# controls aspect ratio of blocks
ocr.trigger = 0.25

# z1 layer and target depth in metres
ocr.z1_layer = 20
ocr.z_target_depth = 10000
ocr.num_z_pad_cells = 10
ocr.z_bottom = 100000
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
