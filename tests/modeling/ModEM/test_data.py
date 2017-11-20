import difflib
import glob
import os
import shutil
import sys
from unittest import TestCase

import matplotlib.pyplot as plt

from examples.create_modem_input import select_periods
from mtpy.modeling.modem import Data
# patch that changes the matplotlib behaviour
from tests import make_temp_dir
from tests.modeling import show_patcher

plt.ion()  # enable interactive
# plt.ioff()  # disable interactive, which will also disable this patch


plt.show = show_patcher(plt.show)


# end of patch

class TestData(TestCase):
    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):
        # for each test, setup a different output dir
        self._output_dir = make_temp_dir(self._testMethodName, base_dir=self._temp_dir)

        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.normpath(
            os.path.join(
                os.path.join(self._temp_dir, 'expected_data_output'),
                self._testMethodName
            )
        )
        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None


edi_paths = [
    "data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "data/edifiles2",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
]
# epsg to project to. Google epsg 'your projection'
epsg_code = 28354
epsg_code = 3112

error_types = [
    # (test_name, error_type_tipper, error_tpye_z, component_error_type)
    ('abs_floor',         'abs',   'floor',         None),
    ('floor_egbert',      'floor', 'egbert',        None),
    ('abs_egbert_floor',  'abs',   'egbert_floor',  None),
    ('floor_mean_od',     'floor', 'mean_od',       None),
    ('abs_mean_od_floor', 'abs',   'mean_od_floor', None),
    ('floor_eigen',       'floor', 'eigen',         None),
    ('abs_eigen_floor',   'abs',   'eigen_floor',   None),
    ('floor_median',      'floor', 'median',        None),
    ('abs_median_floor',  'abs',   'median_floor',  None)
]


def _test_gen(edi_path, error_type_z, error_type_tipper, comp_error_type):
    """
    generate list of tests for the given edi path
    :param index:
    :param edi_path:
    :return:
    """

    def test_func(self):
        if not os.path.isdir(edi_path):
            # input file does not exist, skip test after remove the output dir
            os.rmdir(self._output_dir)
            self.skipTest("edi path does not exist: {}".format(edi_path))

        edi_list = glob.glob(edi_path + '/*.edi')
        period_list = select_periods(edi_list)
        datob = Data(edi_list=edi_list,
                     inv_mode='1',
                     period_list=period_list,
                     epsg=epsg_code,
                     error_type_tipper=error_type_tipper,
                     error_type_z=error_type_z,
                     comp_error_type=comp_error_type,
                     error_floor=10)
        datob.write_data_file(save_path=self._output_dir)

        # check the output
        if self._expected_output_dir:
            output_data_file = os.path.normpath(os.path.join(self._output_dir, "ModEM_Data.dat"))
            self.assertTrue(os.path.isfile(output_data_file), "output data file does not exist")
            expected_data_file = os.path.normpath(os.path.join(self._expected_output_dir,
                                                               "ModEM_Data.dat"))
            self.assertTrue(
                os.path.isfile(expected_data_file),
                "expected output data file does not exist, nothing to compare"
            )

            print "\ncomparing", output_data_file, "and", expected_data_file
            with open(output_data_file, 'r') as output:
                with open(expected_data_file, 'r') as expected:
                    diff = difflib.unified_diff(
                        expected.readlines(),
                        output.readlines(),
                        fromfile='expected',
                        tofile='output'
                    )
                    count = 0
                    for line in diff:
                        sys.stdout.write(line)
                        count += 1
                    self.assertTrue(count == 0, "output different!")
        else:
            print "no expected output exist, nothing to compare"

    return test_func


# generate tests
for edi_path in edi_paths:
    for name, error_type_tipper, error_type_z, comp_error_type in error_types:
        test_func = _test_gen(edi_path, error_type_tipper, error_type_z, comp_error_type)
        test_func.__name__ = "test_{}_{}".format(os.path.basename(edi_path), name)
        setattr(TestData, test_func.__name__, test_func)

if 'test_func' in globals():
    del globals()['test_func']
