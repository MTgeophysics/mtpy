import glob
import os
import shutil
from unittest import TestCase

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

error_types = ['floor']


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
