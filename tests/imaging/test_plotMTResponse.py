import glob
import os
import unittest
from unittest import TestCase

import matplotlib
import matplotlib.pyplot as plt
import shutil

from mtpy.imaging.plot_mt_response import PlotMTResponse
from tests.imaging import reset_matplotlib, ImageTestCase

edi_paths = [
    "",
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests/data/edifiles2"
]


class TestPlotMTResponse(ImageTestCase):
    def test_plot_01(self):
        edi_path = edi_paths[1]
        self._plot(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_plot_02(self):
        edi_path = edi_paths[2]
        self._plot(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03(self):
        edi_path = edi_paths[3]
        self._plot(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    def test_plot_04(self):
        edi_path = edi_paths[4]
        self._plot(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    def test_plot_05(self):
        edi_path = edi_paths[5]
        self._plot(edi_path)

    def _plot(self, edi_path):
        edi_file_list = glob.glob(os.path.join(edi_path, '*.edi'))
        for edi_file in edi_file_list:
            pt_obj = PlotMTResponse(
                fn=edi_file
            )
            pt_obj.plot()
            plt.pause(.5)
