# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:18:12 2017

@author: u64125
"""

from mtpy.modeling.modem import Model
from mtpy.modeling.modem import Data
from mtpy.modeling.modem import Covariance
from mtpy.utils.calculator import get_period_list
from tests import EDI_DATA_DIR2, TEST_TEMP_DIR
import numpy as np

# path to save to
workdir = TEST_TEMP_DIR.joinpath("ModEM")
if not workdir.exists():
    workdir.mkdir()


## period list (won't include periods outside of the range of the edi file) ###
## comment/uncomment your desired method ######################################
###############################################################################

## example to specify start, stop and total number of periods
# start_period = 0.001
# stop_period = 1000
# n_periods = 25
# period_list = np.logspace(np.log10(start_period),
#                          np.log10(stop_period),
#                          n_periods)


# example to specify a number of periods per decade
start_period = 0.002
stop_period = 2000
periods_per_decade = 4
period_list = get_period_list(
    start_period, stop_period, periods_per_decade, include_outside_range=True
)


## an example to use the periods from a particular edi file
# edifile_periods = op.join(edipath,'Synth00.edi')
# eobj = Edi(edifile_periods)
# period_list = 1./eobj.freq

###############################################################################

# list of edi files, search for all files ending with '.edi'
edi_list = list(EDI_DATA_DIR2.glob("*.edi"))


do = Data(
    edi_list=edi_list,
    inv_mode="1",
    save_path=workdir,
    period_list=period_list,
    period_buffer=2,  # factor to stretch interpolation by. For example: if period_buffer=2
    # then interpolated data points will only be included if they are
    # within a factor of 2 of a true data point
    error_type_z=np.array(
        [
            [
                "floor_percent",
                "floor_egbert",
            ],  # error type, options are 'egbert', 'percent', 'mean_od', 'eigen', 'median', 'off_diagonals'
            ["floor_egbert", "percent"],
        ]
    ),  # add floor to apply it as an error floor
    # can supply a 2 x 2 array for each component or a single value
    error_value_z=np.array(
        [[20.0, 5.0], [5.0, 20.0]]  # error floor value in percent
    ),  # can supply a 2 x 2 array for each component or a single value
    error_type_tipper="floor_abs",  # type of error to set in tipper,
    # floor_abs is an absolute value set as a floor
    error_value_tipper=0.03,
    model_epsg=28354  # model epsg, currently set to utm zone 54.
    # See http://spatialreference.org/ to find the epsg code for your projection
)

# Unlike when writing topography from a file, don't modify the
#  elevation of the Data object as we need the station elevations
#  to create surface model.
do.write_data_file()

# create model file
mo = Model(
    station_locations=do.station_locations,
    cell_size_east=8000,
    cell_size_north=8000,
    pad_north=7,  # number of padding cells in each of the north and south directions
    pad_east=7,  # number of east and west padding cells
    pad_z=6,  # number of vertical padding cells
    pad_stretch_v=1.6,  # factor to increase by in padding cells (vertical)
    pad_stretch_h=1.4,  # factor to increase by in padding cells (horizontal)
    pad_num=3,  # number of constant-width cells to add to outside of model before padding cells start
    # this number is currently multiplied by 1.5 internally
    n_air_layers=10,  # number of air layers, set to 0 to incorporate bathymetry only
    res_model=100,  # halfspace resistivity value for reference model
    n_layers=100,  # total number of z layers, including air
    z1_layer=10,  # first layer thickness
    pad_method="stretch",  # method for calculating padding
    z_mesh_method="new",
    z_target_depth=120000,  # depth to bottom of core model (padding after this depth)
)

mo.make_mesh()
mo.write_model_file(save_path=workdir)

# Add topography from EDI data.
# This function takes the same arguments as `add_topography_to_model2`
#  but instead of a DEM file it takes a Data object.
mo.add_topography_from_data(do)
mo.write_model_file(save_path=workdir)

# update data elevations
do.project_stations_on_topography(mo)

# show the mesh
mo.plot_sealevel_resistivity()


co = Covariance()
co.smoothing_east = 0.4
co.smoothing_north = 0.4
co.smoothing_z = 0.4
co.write_covariance_file(model_fn=mo.model_fn)
