# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    From a set of EDI files, create input files for MODEM inversion program.
    With options to modify the resistivity values of horizontal slices at specified depth.

CreationDate:   17/12/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     17/12/2019   FZ initiated this script from the scripts/ModEM_build_inputfiles.py
    LastUpdate:     dd/mm/yyyy
"""

import os
import os.path as op
import sys
import numpy as np
from mtpy.modeling.modem import Model
from mtpy.modeling.modem import Data
from mtpy.modeling.modem import Covariance
from mtpy.core.edi import Edi
from mtpy.utils.calculator import get_period_list


def build_modem_inputfiles(edi_path, output_path):
    """
    For a given se to edi files in edi_path, build the ModEM input files into output_path
    :param edi_path: # path where edi files are located r'C:\mtpywin\mtpy\examples\data\edi_files_2'
    :param output_path: # path to save to  r'C:\test\ModEM'
    :return:
    """

    # not nece os.chdir(r'C:\mtpywin\mtpy')  # change this path to the path where mtpy is installed

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
    period_list = get_period_list(start_period, stop_period, periods_per_decade,
                                  include_outside_range=True)

    ## an example to use the periods from a particular edi file
    # edifile_periods = op.join(edipath,'Synth00.edi')
    # eobj = Edi(edifile_periods)
    # period_list = 1./eobj.freq

    ###############################################################################
    workdir = output_path
    edipath = edi_path
    # list of edi files, search for all files ending with '.edi'
    edi_list = [op.join(edipath, ff) for ff in os.listdir(edipath) if (ff.endswith('.edi'))]

    # make the save path if it doesn't exist
    if not op.exists(workdir):
        os.mkdir(workdir)

    do = Data(edi_list=edi_list,
              inv_mode='1',
              save_path=workdir,
              period_list=period_list,
              period_buffer=2,  # factor to stretch interpolation by. For example: if period_buffer=2
              # then interpolated data points will only be included if they are
              # within a factor of 2 of a true data point
              error_type_z=np.array([['floor_percent', 'floor_egbert'],
                                     # error type, options are 'egbert', 'percent', 'mean_od', 'eigen', 'median', 'off_diagonals'
                                     ['floor_egbert', 'percent']]),  # add floor to apply it as an error floor
              # can supply a 2 x 2 array for each component or a single value
              error_value_z=np.array([[20., 5.],  # error floor value in percent
                                      [5., 20.]]),  # can supply a 2 x 2 array for each component or a single value
              error_type_tipper='floor_abs',  # type of error to set in tipper,
              # floor_abs is an absolute value set as a floor
              error_value_tipper=.03,
              model_epsg=28354  # model epsg, currently set to utm zone 54.
              # See http://spatialreference.org/ to find the epsg code for your projection
              )
    do.write_data_file()

    # set elevations to zero as we need to ensure the stations are on the topography
    do.data_array['elev'] = 0.
    do.write_data_file(fill=False)

    # create model file
    mo = Model(station_locations=do.station_locations,
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
               pad_method='stretch',  # method for calculating padding
               z_mesh_method='new',
               z_target_depth=120000  # depth to bottom of core model (padding after this depth)
               )

    mo.make_mesh()
    mo.write_model_file(save_path=workdir)

    # add topography to res model
    # if the number of air layers is zero - bathymetry only will be added.
    # if the number of air layers is nonzero - topography will be added, discretised into that number of cells
    # mo.add_topography_to_model2(r'C:\mtpywin\mtpy\examples\data\AussieContinent_etopo1.asc')
    mo.add_topography_to_model2(r'C:\Githubz\mtpy\examples\data\AussieContinent_etopo1.asc')
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


# ====================================================================================================
# Entry point of this script
# How to run (Example):
# python examples/cmdline/build_ModEM_inputfiles.py examples/data/edi_files_2 /c/Workspace/edi_files_2
# ====================================================================================================
if __name__ == "__main__":

    USAGE = "\n********************************************\n " \
            "Usage:  python %s edi_dir_path output_dir_path "   \
            "\n********************************************\n " % sys.argv[0]

    if len(sys.argv) < 3:
        print(USAGE)
        sys.exit(1)

    build_modem_inputfiles(sys.argv[1], sys.argv[2])
