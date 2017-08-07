import glob
from unittest import TestCase

import os
import shutil
import matplotlib.pyplot as plt

from examples.create_modem_input import select_periods
from mtpy.modeling.modem_data import Data
from mtpy.modeling.modem_model import Model


class TestModel(TestCase):
    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

    def setUp(self):
        # for each test, setup a different output dir
        self._output_dir = os.path.normpath(os.path.join(self._temp_dir, self._testMethodName))
        if os.path.exists(self._output_dir):
            # clear dir if it already exist
            shutil.rmtree(self._output_dir)
        os.mkdir(self._output_dir)

        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.normpath(
            os.path.join(
                os.path.join(self._temp_dir, 'expected_model_output'),
                self._testMethodName
            )
        )
        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None


edi_paths = [
    "tests\\data\\edifiles",
    "examples\\data\\edi2",
    "examples\\data\\edi_files",
    "tests\\data\\edifiles2",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
]
# epsg to project to. Google epsg 'your projection'
epsg_code = 28354
epsg_code = 3112


def _test_gen(index, edi_path):
    def test_func(self):
        if not os.path.isdir(edi_path):
            # input file does not exist, skip test after remove the output dir
            os.rmdir(self._output_dir)
            self.skipTest("edi path does not exist: {}".format(edi_path))

        # generate data
        edi_list = glob.glob(edi_path + '/*.edi')
        period_list = select_periods(edi_list)
        datob = Data(edi_list=edi_list,
                     inv_mode='1',
                     period_list=period_list,
                     epsg=epsg_code,
                     error_type='egbert',
                     comp_error_type=None,
                     error_floor=10)
        datob.write_data_file(save_path=self._output_dir)

        # create mesh grid model object
        model = Model(Data=datob,
                      epsg=epsg_code,
                      cell_size_east=10000, cell_size_north=10000,  # GA_VIC
                      pad_north=8,  # number of padding cells in each of the north and south directions
                      pad_east=8,  # number of east and west padding cells
                      pad_z=8,  # number of vertical padding cells
                      pad_stretch_v=1.5,  # factor to increase by in padding cells (vertical)
                      pad_stretch_h=1.5,  # factor to increase by in padding cells (horizontal)
                      n_airlayers=0,  # number of air layers 0, 10, 20, depend on topo elev height
                      res_model=100,  # halfspace resistivity value for initial reference model
                      n_layers=50,  # total number of z layers, including air and pad_z
                      z1_layer=50,  # first layer thickness metres, depend
                      z_target_depth=500000
                      )
        model.make_mesh()  # the data file will be re-write in this method. No topo elev file used yet
        model.plot_mesh()
        model.plot_mesh_xy()
        model.plot_mesh_xz()

        # write a model file and initialise a resistivity model
        model.write_model_file(save_path=self._output_dir)

    return test_func


# generate tests
for index, edi_path in enumerate(edi_paths):
    test_func = _test_gen(index, edi_path)
    test_func.__name__ = "test_{}_{}".format(index+1, os.path.basename(edi_path))
    setattr(TestModel, test_func.__name__, test_func)
