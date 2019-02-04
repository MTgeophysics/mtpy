

import glob
import os
import unittest
from unittest import TestCase
import filecmp
import matplotlib
import sys

from tests import make_temp_dir
from tests.imaging import plt_wait

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", file=sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt
    plt.ion()

import numpy as np
from geopandas import GeoDataFrame

from mtpy.core.edi_collection import is_num_in_seq, EdiCollection
from mtpy.core.mt import MT

edi_paths = [
    #"../../data/edifiles",
    #"../../examples/data/edi2",
    "../../examples/data/edi_files",
    #"../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    #"../MT_Datasets/GA_UA_edited_10s-10000s",
    #"../../data/edifiles2"
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
        self.edi_files = glob.glob(os.path.normpath(os.path.abspath(os.path.join(self.edi_path, "*.edi"))))

    @classmethod
    def setUpClass(cls):
        plt.clf()
        cls._temp_dir = make_temp_dir(cls.__name__)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_init(self):
        """
        test if the EdiCollection is initialized correctly
        :return:
        """
        edi_files = glob.glob(os.path.normpath(os.path.abspath(os.path.join(self.edi_path, "*.edi"))))
        self.assertTrue(len(edi_files) == self.edi_collection.num_of_edifiles)
        self.assertTrue(set(edi_files) == set(self.edi_collection.edifiles))
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

    def test_get_tensor_tippers(self):
        mto = self.edi_collection.mt_obj_list[0]
        p = 1./mto.Z.freq[0]
        pt_dict_list_with_interp = self.edi_collection.get_phase_tensor_tippers(p)
        pt_dict_list_no_interp = self.edi_collection.get_phase_tensor_tippers(p, interpolate=False)

        # Asserting parity of results from interpolation and that from without
        assert np.allclose(pt_dict_list_with_interp[0]['phi_min'], mto.pt.phimin[0])
        assert np.allclose(pt_dict_list_no_interp[0]['phi_min'], mto.pt.phimin[0])
    # end func

    def test_plot_stations(self):
        self.edi_collection.plot_stations()
        plt_wait(1)

    def test_display_on_basemap(self):
        self.edi_collection.display_on_basemap()
        plt_wait(1)

    def test_display_on_image(self):
        #self.edi_collection.display_on_image()
        plt_wait(1)

    def test_create_mt_station_gdf(self):
        path = make_temp_dir(self.__class__.__name__ + "_mt_station_gdf", base_dir=self._temp_dir)
        self.edi_collection.create_mt_station_gdf(path)

    def test_create_measurement_csv(self):
        path = make_temp_dir(self.__class__.__name__ + "_measurement_csv", base_dir=self._temp_dir)
        self.edi_collection.create_measurement_csv(path)

    def test_create_phase_tensor_csv(self):
        path = make_temp_dir(self.__class__.__name__ + "_phase_tensor_csv", base_dir=self._temp_dir)
        self.edi_collection.create_phase_tensor_csv(path)

    def test_create_phase_tensor_csv_with_image(self):
        path2 = make_temp_dir(self.__class__.__name__ + "_phase_tensor_csv_with_image", base_dir=self._temp_dir)
        self.edi_collection.create_phase_tensor_csv_with_image(path2)

    def test_export_edi_files(self):
        path3 = make_temp_dir(self.__class__.__name__ + "_export_edi_files_no_interp", base_dir=self._temp_dir)
        path4 = make_temp_dir(self.__class__.__name__ + "_export_edi_files_interp", base_dir=self._temp_dir)

        # Note that this test relies on the fact that the input edi files all have the same periods.
        # EDI files generated based on the periods as in the input edi files should be identical
        # to those produced with interpolation turned on or off.

        # no interp
        self.edi_collection.export_edi_files(path3)

        # interp
        mto = self.edi_collection.mt_obj_list[0]
        plist = 1./mto.Z.freq
        self.edi_collection.export_edi_files(path4, period_list=plist)

        for fn in glob.glob('%s/*.edi'%(path3)):
            f1 = os.path.join(path3, fn)
            f2 = os.path.join(path4, fn)

            assert filecmp.cmp(f1, f2)
        # end for

class _TestFromFile(_BaseTest):
    def setUp(self):
        _BaseTest.setUp(self)
        self.edi_collection = EdiCollection(self.edi_files)


class _TestFromMTObj(_BaseTest):
    def setUp(self):
        _BaseTest.setUp(self)
        mt_objs = [MT(edi_file) for edi_file in self.edi_files]
        self.edi_collection = EdiCollection(mt_objs=mt_objs)


for edi_path in edi_paths:
    if os.path.isdir(edi_path):
        cls_name = "TestEdiCollectionFromFile_%s" % (os.path.basename(edi_path))
        globals()[cls_name] = type(cls_name, (_TestFromFile, unittest.TestCase), {
            "edi_path": edi_path
        })
        cls_name = "TestEdiCollectionFromMTObj_%s" % (os.path.basename(edi_path))
        globals()[cls_name] = type(cls_name, (_TestFromMTObj, unittest.TestCase), {
            "edi_path": edi_path
        })
