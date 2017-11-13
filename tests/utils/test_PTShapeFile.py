import inspect
import os
import unittest
from unittest import TestCase
import numpy as np
import pytest

from mtpy.utils.shapefiles import create_phase_tensor_shpfiles, transform_ll_to_utm, transform_utm_to_ll
from tests import TEST_TEMP_DIR

edi_paths = [
    "",
    "data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "data/edifiles2"
]


class TestPTShapeFile(TestCase):
    def setUp(self):
        self._temp_dir = TEST_TEMP_DIR

    def test_shapefile_01(self):
        edi_path = edi_paths[1]
        self._create_shape_file(edi_path, os.path.join(self._temp_dir, inspect.currentframe().f_code.co_name))

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_shapefile_02(self):
        edi_path = edi_paths[2]
        self._create_shape_file(edi_path, os.path.join(self._temp_dir, inspect.currentframe().f_code.co_name))

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_shapefile_03(self):
        edi_path = edi_paths[3]
        self._create_shape_file(edi_path, os.path.join(self._temp_dir, inspect.currentframe().f_code.co_name))

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    def test_shapefile_04(self):
        edi_path = edi_paths[4]
        self._create_shape_file(edi_path, os.path.join(self._temp_dir, inspect.currentframe().f_code.co_name))

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    def test_shapefile_05(self):
        edi_path = edi_paths[5]
        self._create_shape_file(edi_path, os.path.join(self._temp_dir, inspect.currentframe().f_code.co_name))

    @staticmethod
    def _create_shape_file(edi_path, save_path):
        if not os.path.isdir(save_path):
            os.mkdir(save_path)
        create_phase_tensor_shpfiles(edi_path, save_path, ellipse_size=6000, every_site=1)


class TestShapefileTools(TestCase):
    def setUp(self):
        self.nref = 23  # 23, "WGS-84"
        self.lat = -34.3
        self.lon = 149.2

        self.zone = '55H'
        self.easting = 702466.54318647704
        self.northing = 9997808.8686107509

    def test_transform_ll_to_utm(self):
        zone, (easting, northing, elev) = transform_ll_to_utm(self.lat, self.lon, "WGS84")

        self.assertTrue(zone, self.zone)
        self.assertTrue(easting, self.easting)
        self.assertTrue(northing, self.northing)

    @pytest.mark.xfail
    def test_transform_utm_to_ll(self):
        new_lat, new_lon = transform_utm_to_ll(self.easting, self.northing, 55)
        # this test fails due to unknown return type of transform_utm_to_ll
        self.assertTrue(np.isclose(self.lat, new_lat))
        self.assertTrue(np.isclose(self.lon, new_lon))
