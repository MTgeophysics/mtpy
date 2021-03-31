from unittest import TestCase
import numpy as np
import os

from mtpy.modeling.modem import Model, Data
from mtpy.utils.calculator import make_log_increasing_array
from tests import make_temp_dir, SAMPLE_DIR
from tests.modeling import diff_files


class TestModEM_Model(TestCase):
    def setUp(self):
        self.model_epsg = 28355
        # self._temp_dir = make_temp_dir(self.__name__)
        self._temp_dir = make_temp_dir("TestModEM_Model_tmp")
        # directory to save created input files
        self._output_dir = make_temp_dir(self._testMethodName, base_dir=self._temp_dir)
        self._model_dir = os.path.join(SAMPLE_DIR, "ModEM")
        self._sgrid_fn = os.path.join(SAMPLE_DIR, "gocad", "ModEM_Model_File.sg")
        self._model_fn = os.path.join(self._model_dir, "ModEM_Model_File.rho")
        self._model_fn_old_z_mesh = os.path.join(
            SAMPLE_DIR, "ModEM_old_make_z_mesh", "ModEM_Model_File.rho"
        )
        self._data_fn = os.path.join(self._model_dir, "ModEM_Data.dat")

    def test_read_gocad_sgrid_file(self):

        if not os.path.isdir(self._model_dir):
            self._model_dir = None

        output_fn = os.path.basename(self._model_fn)

        # read data file to get centre position
        dObj = Data()
        dObj.read_data_file(data_fn=self._data_fn)

        # create a model object using the data object and read in gocad sgrid file
        mObj = Model(data_obj=dObj, save_path=self._output_dir)
        mObj.read_gocad_sgrid_file(self._sgrid_fn)
        mObj.write_model_file()

        output_data_file = os.path.normpath(os.path.join(self._output_dir, output_fn))

        self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

        expected_data_file = os.path.normpath(self._model_fn_old_z_mesh)

        self.assertTrue(
            os.path.isfile(expected_data_file),
            "Ref output data file does not exist, nothing to compare with",
        )

        is_identical, msg = diff_files(output_data_file, expected_data_file)
        print(msg)
        self.assertTrue(
            is_identical, "The output file is not the same with the baseline file."
        )

    def test_write_gocad_sgrid_file(self):

        if not os.path.exists(self._sgrid_fn):
            self._sgrid_fn = None

        output_fn = os.path.basename(self._sgrid_fn)

        # read data file to get centre position
        dObj = Data()
        dObj.read_data_file(data_fn=self._data_fn)

        # get centre coordinates
        centre = np.array([0.0, 0.0, 0.0])
        centre[0] = dObj.center_point["east"]
        centre[1] = dObj.center_point["north"]

        # create a model object using the data object and read in gocad sgrid file
        mObj = Model(data_obj=dObj)
        mObj.read_model_file(model_fn=self._model_fn)
        mObj.save_path = self._output_dir
        mObj.write_gocad_sgrid_file(
            origin=centre, fn=os.path.join(self._output_dir, output_fn[:-3])
        )

        output_data_file = os.path.normpath(os.path.join(self._output_dir, output_fn))

        self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

        expected_data_file = os.path.normpath(self._sgrid_fn)

        self.assertTrue(
            os.path.isfile(expected_data_file),
            "Ref output data file does not exist, nothing to compare with",
        )

        is_identical, msg = diff_files(output_data_file, expected_data_file)
        print(msg)
        self.assertTrue(
            is_identical, "The output file is not the same with the baseline file."
        )

    def test_make_z_mesh_new(self):

        z1_layer = 10
        z_target_depth = 5000
        n_layers = 30
        n_airlayers = 10
        pad_z = 4
        pad_stretch_v = 1.4

        mObj = Model(
            z1_layer=z1_layer,
            z_target_depth=z_target_depth,
            n_layers=n_layers,
            n_air_layers=n_airlayers,
            pad_z=pad_z,
            pad_stretch_v=pad_stretch_v,
        )
        z_nodes, z_grid = mObj.make_z_mesh_new()

        expected_air_layers = [
            10.0,
            10.0,
            10.0,
            20.0,
            20.0,
            20.0,
            30.0,
            30.0,
            40.0,
            50.0,
        ]
        self.assertTrue(np.all(z_nodes[:n_airlayers] == expected_air_layers))
        expected_air_layers = [
            0.0,
            10.0,
            20.0,
            30.0,
            50.0,
            70.0,
            90.0,
            120.0,
            150.0,
            190.0,
            240.0,
        ]
        self.assertTrue(np.all(z_grid[: n_airlayers + 1] == expected_air_layers))

        # check core model part
        testnodes10_5000_16 = [
            60.0,
            70.0,
            80.0,
            100.0,
            100.0,
            100.0,
            200.0,
            200.0,
            200.0,
            300.0,
            300.0,
            400.0,
            500.0,
            600.0,
            700.0,
            800.0,
        ]
        testgrid10_5000_16 = [
            240.0,
            300.0,
            370.0,
            450.0,
            550.0,
            650.0,
            750.0,
            950.0,
            1150.0,
            1350.0,
            1650.0,
            1950.0,
            2350.0,
            2850.0,
            3450.0,
            4150.0,
            4950.0,
        ]
        self.assertTrue(np.all(z_nodes[n_airlayers:-pad_z] == testnodes10_5000_16))
        self.assertTrue(np.all(z_grid[n_airlayers:-pad_z] == testgrid10_5000_16))

        # check padding part
        testnodespad = np.around(
            testnodes10_5000_16[-1] * (pad_stretch_v ** np.arange(1, pad_z + 1)), -2
        )
        testgridpad = (
            np.array([testnodespad[:i].sum() for i in range(1, pad_z + 1)])
            + testgrid10_5000_16[-1]
        )

        self.assertTrue(np.all(z_nodes[-pad_z:] == testnodespad))
        self.assertTrue(np.all(z_grid[-pad_z:] == testgridpad))
