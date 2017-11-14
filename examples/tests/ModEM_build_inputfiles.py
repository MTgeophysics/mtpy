# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:18:12 2017

@author: u64125
"""
import os
os.chdir(r'C:\Git\mtpy') # change this path to the path where mtpy is installed
import os.path as op
from mtpy.modeling.ModEM import Data, Model
from mtpy.modeling.modem import Covariance
import mtpy.core.edi as mtedi
import numpy as np
import matplotlib.pyplot as plt

# path to save to
workdir = r'C:\Git\mtpy\examples\model_files\ModEM'

# path where edi files are located
edipath = r'C:\Git\mtpy\examples\data\edi_files'

# period list (will not include periods outside of the range of the edi file)
start_period = -2
stop_period = 3
n_periods = 17
period_list = np.logspace(start_period,stop_period,n_periods)


# list of edi files, search for all files ending with '.edi'
edi_list = [op.join(edipath,ff) for ff in os.listdir(edipath) if (ff.endswith('.edi'))]

if not op.exists(workdir):
    os.mkdir(workdir)


do = Data(edi_list=edi_list,
               inv_mode = '1',
               save_path=workdir,
               period_list=period_list,
               error_type_z='floor_egbert',
               error_value_z=5,
               error_type_tipper = 'floor_abs',
               error_value_tipper =.03,
               model_epsg=28354 # model epsg, currently set to utm zone 54
               )
do.write_data_file()
do.data_array['elev'] = 0.
do.write_data_file(fill=False)

# create model file
mo = Model(station_locations=do.station_locations,
                cell_size_east=500,
                cell_size_north=500,
                pad_north=7, # number of padding cells in each of the north and south directions
                pad_east=7,# number of east and west padding cells
                pad_z=6, # number of vertical padding cells
                pad_stretch_v=1.6, # factor to increase by in padding cells (vertical)
                pad_stretch_h=1.4, # factor to increase by in padding cells (horizontal)
                n_airlayers = 10, #number of air layers
                res_model=100, # halfspace resistivity value for reference model
                n_layers=90, # total number of z layers, including air
                z1_layer=10, # first layer thickness
                pad_method='stretch',
                z_target_depth=120000)

mo.make_mesh()
mo.write_model_file(save_path=workdir)

# add topography to res model
mo.add_topography_to_model2(r'C:\Git\mtpy\examples\data\AussieContinent_etopo1.asc')
mo.write_model_file(save_path=workdir)

do.project_stations_on_topography(mo)


co = Covariance()
co.write_covariance_file(model_fn=mo.model_fn)
