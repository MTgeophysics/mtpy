# -*- coding: utf-8 -*-
"""
This script is to create input files required to run occam1d inversion, which include:

OccamStartup1D
Model1D
Occam1d_DataFile_DET.dat
These input files are created from a standard edi data file.

"""

import os.path as op
import os
import mtpy.modeling.occam1d as mtoc1d  # Wrapper class to interact with Occam1D

# directory to save created input files
savepath = r"C:/tmp"
if not op.exists(savepath):
    os.mkdir(savepath)

# edi path and file name
edi_file = r"C:/mtpywin/mtpy/examples/data/edi_files_2/Synth00.edi"

# create data file
ocd = mtoc1d.Data()  # create an object and assign values to arguments

ocd.write_data_file(
    edi_file=edi_file,
    mode="te",  # det mode
    save_path=savepath,
    res_errorfloor=4,  # error floor in percentage
    phase_errorfloor=2,  # error floor in degrees
    remove_outofquadrant=True,
)

# create model file
ocm = mtoc1d.Model(
    n_layers=100,  # number of layers
    target_depth=40000,  # target depth in metres, before padding
    bottom_layer=100000,
    z1_layer=10,  # first layer thickness in metres
)
ocm.write_model_file(save_path=savepath)

# create startup file
ocs = mtoc1d.Startup(
    data_fn=ocd.data_fn,  # basename of data file *default* is Occam1DDataFile
    model_fn=ocm.model_fn,  # basename for model file *default* is Model1D
    max_iter=200,  # maximum number of iterations to run
    target_rms=1.0,
)
ocs.write_startup_file()
