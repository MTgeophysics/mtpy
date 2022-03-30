# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    testing ModEM input files builder modules.

    this test suite create all output for modem and compare the created files with the baseline files

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
from pathlib import Path
import os
from unittest import TestCase

import numpy as np

from mtpy.modeling.modem import Data, Model, Covariance
from tests import EDI_DATA_DIR, EDI_DATA_DIR2, AUS_TOPO_FILE, SAMPLE_DIR, make_temp_dir
from tests.modeling import diff_files
from mtpy.utils.calculator import get_period_list


class TestModemInputFilesBuilder(TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):

        # directory to save created input files
        self._output_dir = make_temp_dir(self._testMethodName, base_dir=self._temp_dir)

        self._expected_output_dir = Path(SAMPLE_DIR, "ModEM")
        if not self._expected_output_dir.is_dir():
            self._expected_output_dir = None

    def test_fun(self):

        edipath = EDI_DATA_DIR  # path where edi files are located
        # set the dir to the output from the previously correct run
        self._expected_output_dir = Path(SAMPLE_DIR, "ModEM")

        # period list (will not include periods outside of the range of the edi file)
        start_period = -2
        stop_period = 3
        n_periods = 17
        period_list = np.logspace(start_period, stop_period, n_periods)

        # list of edi files, search for all files ending with '.edi'
        edi_list = list(edipath.glob("*.edi"))

        do = Data(
            edi_list=edi_list,
            inv_mode="1",
            save_path=self._output_dir,
            period_list=period_list,
            error_type_z="floor_egbert",
            error_value_z=5,
            error_type_tipper="floor_abs",
            error_value_tipper=0.03,
            model_epsg=28354,  # model epsg, currently set to utm zone 54
        )

        # BM: write here to fill the data object, but this is not the data file being compared.
        do.write_data_file()

        # create model file
        mo = Model(
            station_locations=do.station_locations,
            cell_size_east=500,
            cell_size_north=500,
            pad_north=7,  # number of padding cells in each of the north and south directions
            pad_east=7,  # number of east and west padding cells
            pad_z=6,  # number of vertical padding cells
            pad_stretch_v=1.6,  # factor to increase by in padding cells (vertical)
            pad_stretch_h=1.4,  # factor to increase by in padding cells (horizontal)
            n_air_layers=10,  # number of air layers
            res_model=100,  # halfspace resistivity value for reference model
            n_layers=90,  # total number of z layers, including air
            z1_layer=10,  # first layer thickness
            pad_method="stretch",
            z_target_depth=120000,
        )

        mo.make_mesh()
        mo.write_model_file(save_path=self._output_dir)

        # add topography to res model
        mo.add_topography_to_model2(AUS_TOPO_FILE)
        # mo.add_topography_to_model2(r'E:/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc')
        mo.write_model_file(save_path=self._output_dir)

        # BM: note this function makes a call to `write_data_file` and this is the datafile being
        #  compared!
        do.project_stations_on_topography(mo)

        co = Covariance()
        co.write_covariance_file(model_fn=mo.model_fn)

        # BM: if this test is failing check that the correct filenames are being selected
        #   for comparison
        for test_output, expected_output in (
            ("ModEM_Data_topo.dat", "ModEM_Data.dat"),
            ("covariance.cov", "covariance.cov"),
            ("ModEM_Model_File.rho", "ModEM_Model_File.rho"),
        ):
            output_data_file = Path(self._output_dir, test_output)

            self.assertTrue(output_data_file.is_file(), "output data file not found")

            expected_data_file = Path(self._expected_output_dir, expected_output)

            self.assertTrue(
                expected_data_file.is_file(),
                "Ref output data file does not exist, nothing to compare with",
            )

            # print ("Comparing", output_data_file, "and", expected_data_file)

            is_identical, msg = diff_files(output_data_file, expected_data_file)
            print(msg)
            self.assertTrue(
                is_identical, "The output file is not the same with the baseline file."
            )

    def test_fun_edi_elevation(self):

        edipath = EDI_DATA_DIR  # path where edi files are located
        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.join(SAMPLE_DIR, "ModEM")

        # period list (will not include periods outside of the range of the edi file)
        start_period = -2
        stop_period = 3
        n_periods = 17
        period_list = np.logspace(start_period, stop_period, n_periods)

        # list of edi files, search for all files ending with '.edi'
        edi_list = [
            os.path.join(edipath, ff)
            for ff in os.listdir(edipath)
            if (ff.endswith(".edi"))
        ]

        do = Data(
            edi_list=edi_list,
            inv_mode="1",
            save_path=self._output_dir,
            period_list=period_list,
            error_type_z="floor_egbert",
            error_value_z=5,
            error_type_tipper="floor_abs",
            error_value_tipper=0.03,
            model_epsg=28354,  # model epsg, currently set to utm zone 54
        )

        do.write_data_file()

        # create model file
        mo = Model(
            station_locations=do.station_locations,
            cell_size_east=500,
            cell_size_north=500,
            pad_north=7,  # number of padding cells in each of the north and south directions
            pad_east=7,  # number of east and west padding cells
            pad_z=6,  # number of vertical padding cells
            pad_stretch_v=1.6,  # factor to increase by in padding cells (vertical)
            pad_stretch_h=1.4,  # factor to increase by in padding cells (horizontal)
            n_air_layers=10,  # number of air layers
            res_model=100,  # halfspace resistivity value for reference model
            n_layers=90,  # total number of z layers, including air
            z1_layer=10,  # first layer thickness
            pad_method="stretch",
            z_target_depth=120000,
        )

        mo.make_mesh()
        mo.write_model_file(save_path=self._output_dir)

        # Add topography from EDI data
        mo.add_topography_from_data(do)
        mo.write_model_file(save_path=self._output_dir)

        # BM: note this function makes a call to `write_data_file` and this is the datafile being
        #  compared!
        do.project_stations_on_topography(mo)

        co = Covariance()
        co.write_covariance_file(model_fn=mo.model_fn)

        # BM: if this test is failing check that the correct filenames are being selected
        #   for comparison
        for test_output, expected_output in (
            ("ModEM_Data_topo.dat", "ModEM_Data_EDI_elev.dat"),
            ("covariance.cov", "covariance_EDI_elev.cov"),
            ("ModEM_Model_File.rho", "ModEM_Model_File_EDI_elev.rho"),
        ):
            output_data_file = os.path.normpath(
                os.path.join(self._output_dir, test_output)
            )

            self.assertTrue(
                os.path.isfile(output_data_file), "output data file not found"
            )

            expected_data_file = os.path.normpath(
                os.path.join(self._expected_output_dir, expected_output)
            )

            self.assertTrue(
                os.path.isfile(expected_data_file),
                "Ref output data file '{}' does not exist, nothing to compare with".format(
                    expected_data_file
                ),
            )

            is_identical, msg = diff_files(output_data_file, expected_data_file)
            print(msg)
            self.assertTrue(
                is_identical,
                "The output file '{}' is not the same with the baseline file '{}'.".format(
                    output_data_file, expected_data_file
                ),
            )

    def test_fun_rotate(self):
        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.join(SAMPLE_DIR, "ModEM_rotate40")

        edipath = EDI_DATA_DIR2

        # example to specify a number of periods per decade
        start_period = 0.002
        stop_period = 2000
        periods_per_decade = 4
        period_list = get_period_list(
            start_period, stop_period, periods_per_decade, include_outside_range=True
        )

        # list of edi files, search for all files ending with '.edi'
        edi_list = [
            os.path.join(edipath, ff)
            for ff in os.listdir(edipath)
            if (ff.endswith(".edi"))
        ]

        do = Data(
            edi_list=edi_list,
            inv_mode="1",
            save_path=self._output_dir,
            period_list=period_list,
            period_buffer=2,  # factor to stretch interpolation by. For example: if period_buffer=2
            # then interpolated data points will only be included if they are
            # within a factor of 2 of a true data point
            error_type_z="floor_egbert",  # error type (egbert is % of sqrt(zxy*zyx))
            # floor means apply it as an error floor
            error_value_z=5,  # error floor (or value) in percent
            error_type_tipper="floor_abs",  # type of error to set in tipper,
            # floor_abs is an absolute value set as a floor
            error_value_tipper=0.03,
            rotation_angle=40,
            model_epsg=28354  # model epsg, currently set to utm zone 54.
            # See http://spatialreference.org/ to find the epsg code for your projection
        )
        do.write_data_file()
        do.data_array["elev"] = 0.0
        do.write_data_file(fill=False)

        # mesh rotation angle is the opposite direction to the rotation of the stations
        if do.rotation_angle == 0:
            mesh_rotation_angle = 0
        else:
            mesh_rotation_angle = -do.rotation_angle

        # create model file
        mo = Model(
            stations_object=do.station_locations,
            cell_size_east=8000,
            cell_size_north=8000,
            pad_north=7,  # number of padding cells in each of the north and south directions
            pad_east=7,  # number of east and west padding cells
            pad_z=6,  # number of vertical padding cells
            pad_stretch_v=1.6,  # factor to increase by in padding cells (vertical)
            pad_stretch_h=1.4,  # factor to increase by in padding cells (horizontal)
            n_air_layers=10,  # number of air layers
            res_model=100,  # halfspace resistivity value for reference model
            n_layers=100,  # total number of z layers, including air
            z1_layer=10,  # first layer thickness
            pad_method="stretch",  # method for calculating padding
            z_mesh_method="new",
            z_target_depth=120000,  # depth to bottom of core model (padding after this depth)
            mesh_rotation_angle=mesh_rotation_angle,
        )

        mo.make_mesh()
        mo.write_model_file(save_path=self._output_dir)

        # add topography to res model
        mo.add_topography_to_model2(AUS_TOPO_FILE)
        mo.write_model_file(save_path=self._output_dir)

        co = Covariance()
        co.smoothing_east = 0.4
        co.smoothing_north = 0.4
        co.smoothing_z = 0.4
        co.write_covariance_file(model_fn=mo.model_fn)

        for afile in ("ModEM_Data.dat", "covariance.cov", "ModEM_Model_File.rho"):
            output_data_file = os.path.normpath(os.path.join(self._output_dir, afile))

            self.assertTrue(
                os.path.isfile(output_data_file), "output data file not found"
            )

            expected_data_file = os.path.normpath(
                os.path.join(self._expected_output_dir, afile)
            )

            self.assertTrue(
                os.path.isfile(expected_data_file),
                "Ref output data file does not exist, nothing to compare with",
            )

            # print ("Comparing", output_data_file, "and", expected_data_file)

            is_identical, msg = diff_files(output_data_file, expected_data_file)
            print(msg)
            self.assertTrue(
                is_identical, "The output file is not the same with the baseline file."
            )
