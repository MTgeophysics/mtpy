import difflib
import glob
import os
import shutil
from unittest import TestCase

import sys

from examples.create_modem_input import select_periods
from mtpy.modeling.modem_data import Data


class TestData(TestCase):
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
                os.path.join(self._temp_dir, 'expected_output'),
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

error_types = ['floor', 'value', 'egbert', 'floor_egbert']


def _test_gen(index, edi_path):
    """
    generate list of tests for the given edi path
    :param index:
    :param edi_path:
    :return:
    """
    tests = []

    for error_type in error_types:
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
                         error_type=error_type,
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

        tests.append(
            (
                error_type,  # test name
                test_func  # test func
            )
        )
    return tests


# generate tests
for index, edi_path in enumerate(edi_paths):
    tests = _test_gen(index, edi_path)
    for name, test_func in tests:
        test_func.__name__ = "test_{}_{}_{}".format(index + 1, os.path.basename(edi_path), name)
        setattr(TestData, test_func.__name__, test_func)
