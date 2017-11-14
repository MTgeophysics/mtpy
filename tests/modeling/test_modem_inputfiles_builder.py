# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    testing ModEM input files builder modules.

References:
    examples/tests/ModEM_build_inputfiles.py

CreationDate:   1/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     07/11/2017   YG
                    01/11/2017   FZ


current error:
    failed then comparing the output files, the last debit at 1e-6 sometimes are different.
"""

# import section

import os

from mtpy.modeling.ModEM import Data, Model
from mtpy.modeling.modem import Covariance
import numpy as np

import shutil
from unittest import TestCase

from tests import TEST_TEMP_DIR, EDI_DATA_DIR, AUS_TOPO_FILE, SAMPLE_DIR
from tests.modeling import _diff_files


class TestModemInputFilesBuilder(TestCase):
    def setUp(self):

        # directory to save created input files
        self._output_dir = os.path.join(TEST_TEMP_DIR, 'ModEM')
        if os.path.exists(self._output_dir):
            # clear dir if it already exist
            shutil.rmtree(self._output_dir)

        os.mkdir(self._output_dir)

        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.join(SAMPLE_DIR, 'ModEM')

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None

    def test_fun(self):

        edipath = EDI_DATA_DIR  # path where edi files are located

        # period list (will not include periods outside of the range of the edi file)
        start_period = -2
        stop_period = 3
        n_periods = 17
        period_list = np.logspace(start_period, stop_period, n_periods)

        # list of edi files, search for all files ending with '.edi'
        edi_list = [os.path.join(edipath, ff) for ff in os.listdir(edipath) if (ff.endswith('.edi'))]

        do = Data(edi_list=edi_list,
                  inv_mode='1',
                  save_path=self._output_dir,
                  period_list=period_list,
                  error_type_z='floor_egbert',
                  error_value_z=5,
                  error_type_tipper='floor_abs',
                  error_value_tipper=.03,
                  model_epsg=28354  # model epsg, currently set to utm zone 54
                  )
        do.write_data_file()
        do.data_array['elev'] = 0.
        do.write_data_file(fill=False)

        # create model file
        mo = Model(station_locations=do.station_locations,
                   cell_size_east=500,
                   cell_size_north=500,
                   pad_north=7,  # number of padding cells in each of the north and south directions
                   pad_east=7,  # number of east and west padding cells
                   pad_z=6,  # number of vertical padding cells
                   pad_stretch_v=1.6,  # factor to increase by in padding cells (vertical)
                   pad_stretch_h=1.4,  # factor to increase by in padding cells (horizontal)
                   n_airlayers=10,  # number of air layers
                   res_model=100,  # halfspace resistivity value for reference model
                   n_layers=90,  # total number of z layers, including air
                   z1_layer=10,  # first layer thickness
                   pad_method='stretch',
                   z_target_depth=120000)

        mo.make_mesh()
        mo.write_model_file(save_path=self._output_dir)

        # add topography to res model
        mo.add_topography_to_model2(AUS_TOPO_FILE)
        # mo.add_topography_to_model2(r'E:/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc')
        mo.write_model_file(save_path=self._output_dir)

        do.project_stations_on_topography(mo)

        co = Covariance()
        co.write_covariance_file(model_fn=mo.model_fn)

        for afile in ("ModEM_Data.dat", "covariance.cov", "ModEM_Model_File.rho"):
            output_data_file = os.path.normpath(os.path.join(self._output_dir, afile))

            self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

            expected_data_file = os.path.normpath(os.path.join(self._expected_output_dir, afile))

            self.assertTrue(os.path.isfile(expected_data_file),
                            "Ref output data file does not exist, nothing to compare with"
                            )

            # print ("Comparing", output_data_file, "and", expected_data_file)

            is_identical, msg = _diff_files(output_data_file, expected_data_file)
            print msg
            self.assertTrue(is_identical, "The output file is not the same with the baseline file.")
