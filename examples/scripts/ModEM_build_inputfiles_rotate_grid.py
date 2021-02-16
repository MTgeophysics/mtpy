# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:18:12 2017

@author: u64125
"""

from mtpy.modeling.modem import Model
from mtpy.modeling.modem import Data
from mtpy.modeling.modem import Covariance
from mtpy.utils.calculator import get_period_list
from tests import EDI_DATA_DIR2, TEST_TEMP_DIR, AUS_TOPO_FILE

# path to save to
workdir = TEST_TEMP_DIR.joinpath("ModEM")
if not workdir.exists():
    workdir.mkdir()

## period list (won't include periods outside of the range of the edi file) ###
## comment/uncomment your desired method ######################################
###############################################################################

## example to specify start, stop and total number of periods
#start_period = 0.001
#stop_period = 1000
#n_periods = 25
#period_list = np.logspace(np.log10(start_period),
#                          np.log10(stop_period),
#                          n_periods)


# example to specify a number of periods per decade
start_period = 0.002
stop_period = 2000
periods_per_decade = 4
period_list = get_period_list(start_period,stop_period,periods_per_decade,
                                 include_outside_range=True)


## an example to use the periods from a particular edi file
#edifile_periods = op.join(edipath,'Synth00.edi')
#eobj = Edi(edifile_periods)
#period_list = 1./eobj.freq

###############################################################################

# list of edi files, search for all files ending with '.edi'
edi_list = list(EDI_DATA_DIR2.glob("*.edi"))

do = Data(edi_list=edi_list,
               inv_mode = '1',
               save_path=workdir,
               period_list=period_list,
               period_buffer = 2, # factor to stretch interpolation by. For example: if period_buffer=2
                                 # then interpolated data points will only be included if they are
                                 # within a factor of 2 of a true data point
               error_type_z='floor_egbert', # error type (egbert is % of sqrt(zxy*zyx))
                                            # floor means apply it as an error floor
               error_value_z=5, # error floor (or value) in percent
               error_type_tipper = 'floor_abs', # type of error to set in tipper, 
                                                # floor_abs is an absolute value set as a floor
               error_value_tipper =.03,
               rotation_angle = -45,
               model_epsg=28354 # model epsg, currently set to utm zone 54. 
                                # See http://spatialreference.org/ to find the epsg code for your projection
               )
do.write_data_file()
do.data_array['elev'] = 0.
do.write_data_file(fill=False)

# create model file
mo = Model(stations_object=do.station_locations,
                cell_size_east=8000,
                cell_size_north=8000,
                pad_north=7, # number of padding cells in each of the north and south directions
                pad_east=7,# number of east and west padding cells
                pad_z=6, # number of vertical padding cells
                pad_stretch_v=1.6, # factor to increase by in padding cells (vertical)
                pad_stretch_h=1.4, # factor to increase by in padding cells (horizontal)
                n_air_layers = 10, #number of air layers
                res_model=100, # halfspace resistivity value for reference model
                n_layers=100, # total number of z layers, including air
                z1_layer=10, # first layer thickness
                pad_method='stretch', # method for calculating padding
                z_mesh_method='new',
                z_target_depth=120000, # depth to bottom of core model (padding after this depth)
                )

mo.make_mesh()
mo.write_model_file(save_path=workdir)
mo.plot_mesh()

## add topography to res model
mo.add_topography_to_model2(AUS_TOPO_FILE)
mo.plot_topography()
mo.write_model_file(save_path=workdir)

# update data elevations
do.project_stations_on_topography(mo)


co = Covariance()
co.smoothing_east = 0.4
co.smoothing_north = 0.4
co.smoothing_z = 0.4
co.write_covariance_file(model_fn=mo.model_fn)
