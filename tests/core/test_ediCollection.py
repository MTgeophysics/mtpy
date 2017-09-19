import glob
import os
import unittest
from unittest import TestCase

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive Agg backend")
    # matplotlib.use('Agg')
    test_image = False
    import matplotlib.pyplot as plt
else:
    test_image = True
    # matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.ion()

import numpy as np
from geopandas import GeoDataFrame

from mtpy.core.edi_collection import is_num_in_seq, EdiCollection
from mtpy.core.mt import MT

edi_paths = [
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests/data/edifiles2"
]


class TestUtilities(TestCase):
    def test_is_num_in_seq(self):
        self.assertTrue(is_num_in_seq(3, [1., 2., 3.]))
        self.assertTrue(is_num_in_seq(2, [1., 2.19999999, 3], atol=.2))
        self.assertTrue(is_num_in_seq(2, [1.811111, 4, 3], atol=.2))
        self.assertFalse(is_num_in_seq(2, [1, 3, 4]))
        # self.assertFalse(is_num_in_seq(2, [1, 2.2, 3.2], atol=.2))
        # self.assertFalse(is_num_in_seq(2, [1.8, 4, 3.2], atol=.2))
        # check boundary condition
        self.assertTrue(is_num_in_seq(1, [0, 1.09999999, 2], atol=0.1))
        self.assertTrue(is_num_in_seq(1, [0, 0.90000001, 2], atol=0.1))
        # the next two assertions may be different from implementation
        # due to how python represents floating point numbers
        self.assertFalse(is_num_in_seq(1, [0, 1.1, 2], atol=0.1))
        self.assertTrue(is_num_in_seq(1, [0, 0.9, 2], atol=0.1))
        self.assertFalse(is_num_in_seq(1, [0, 1.10000001, 2], atol=.1))
        self.assertFalse(is_num_in_seq(1, [0, 0.89999999, 2], atol=.1))


class _BaseTest(object):
    def setUp(self):
        self.edi_files = glob.glob(os.path.join(self.edi_path, "*.edi"))

    @classmethod
    def setUpClass(cls):
        plt.clf()
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_init(self):
        """
        test if the EdiCollection is initialized correctly
        :return:
        """
        edi_files = glob.glob(os.path.join(self.edi_path, "*.edi"))
        self.assertTrue(edi_files == self.edi_collection.edifiles)
        self.assertTrue(len(edi_files) == self.edi_collection.num_of_edifiles)
        self.assertTrue(self.edi_collection.ptol == 0.05)  # default
        self.assertTrue(isinstance(self.edi_collection.all_frequencies, list))  # not none and none-empty
        self.assertTrue(isinstance(self.edi_collection.mt_periods, np.ndarray))
        self.assertTrue(isinstance(self.edi_collection.all_unique_periods, np.ndarray))
        self.assertTrue(isinstance(self.edi_collection.geopdf, GeoDataFrame))
        self.assertTrue(isinstance(self.edi_collection.bound_box_dict, dict))

    def test_get_periods_by_stats(self):
        percentages = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0]
        periods = []
        for percent in percentages:
            new_periods = set(self.edi_collection.get_periods_by_stats(percentage=percent))
            for sel_periods in periods:
                if not sel_periods.issubset(new_periods):
                    self.fail()

            periods.append(new_periods)

    def test_plot_stations(self):
        if test_image:
            self.edi_collection.plot_stations()
            plt.pause(1)
        else:
            self.skipTest("Skipped due to matplotlib issues")

    def test_display_on_basemap(self):
        if test_image:
            self.edi_collection.display_on_basemap()
            plt.pause(1)
        else:
            self.skipTest("Skipped due to matplotlib issues")

    def test_display_on_image(self):
        if test_image:
            self.edi_collection.display_on_image()
            plt.pause(1)
        else:
            self.skipTest("Skipped due to matplotlib issues")

    def test_create_mt_station_gdf(self):
        path = os.path.join(self._temp_dir, self.__class__.__name__ + "_mt_station_gdf")
        if not os.path.exists(path):
            os.mkdir(path)
        self.edi_collection.create_mt_station_gdf(path)

    def test_create_measurement_csv(self):
        path = os.path.join(self._temp_dir, self.__class__.__name__ + "_measurement_csv")
        if not os.path.exists(path):
            os.mkdir(path)
        self.edi_collection.create_measurement_csv(path)

    def test_create_phase_tensor_csv(self):
        path = os.path.join(self._temp_dir, self.__class__.__name__ + "_phase_tensor_csv")
        if not os.path.exists(path):
            os.mkdir(path)
        self.edi_collection.create_phase_tensor_csv(path)

    def test_create_phase_tensor_csv_with_image(self):
        if test_image:
            path2 = os.path.join(self._temp_dir, self.__class__.__name__ + "_phase_tensor_csv_with_image")
            if not os.path.exists(path2):
                os.mkdir(path2)
            self.edi_collection.create_phase_tensor_csv_with_image(path2)
        else:
            self.skipTest("Skipped due to matplotlib issues")


class TsetFromFile(_BaseTest):
    def setUp(self):
        _BaseTest.setUp(self)
        self.edi_collection = EdiCollection(self.edi_files)


class TestFromMTObj(_BaseTest):
    def setUp(self):
        _BaseTest.setUp(self)
        mt_objs = [MT(edi_file) for edi_file in self.edi_files]
        self.edi_collection = EdiCollection(mt_objs=mt_objs)


for edi_path in edi_paths:
    if os.path.isdir(edi_path):
        cls_name = "TestEdiCollectionFromFile_%s" % (os.path.basename(edi_path))
        globals()[cls_name] = type(cls_name, (TsetFromFile, unittest.TestCase), {
            "edi_path": edi_path
        })
        cls_name = "TestEdiCollectionFromMTObj_%s" % (os.path.basename(edi_path))
        globals()[cls_name] = type(cls_name, (TestFromMTObj, unittest.TestCase), {
            "edi_path": edi_path
        })

if 'TsetFromFile' in globals():
    del globals()['TsetFromFile']
if 'TestFromMTObj' in globals():
    del globals()['TestFromMTObj']
