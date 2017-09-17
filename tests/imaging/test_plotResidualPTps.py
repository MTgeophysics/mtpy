import glob
import inspect
import unittest
from unittest import TestCase

# configure matplotlib for testing
import matplotlib.pyplot as plt
import os

from mtpy.imaging.plotresidualptps import PlotResidualPTps

plt.ion()

edi_paths = [
    "",
    "tests\\data\\edifiles",
    "examples\\data\\edi2",
    "examples\\data\\edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests\\data\\edifiles2"
]


class TestPlotResidualPTps(TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def setUp(self):
        plt.clf()

    def test_plot_01(self):
        edi_path_1 = edi_paths[1]
        edi_path_2 = edi_paths[2]
        self._plot(edi_path_1, edi_path_2, "%s.png" % inspect.currentframe().f_code.co_name)

    def _plot(self, edi_path_1, edi_path_2, save_figure_path):
        edi_file_list_1 = glob.glob(os.path.join(edi_path_1, "*.edi"))
        edi_file_list_2 = glob.glob(os.path.join(edi_path_2, "*,edi"))
        save_figure_path = os.path.join(self._temp_dir, save_figure_path)
        plt_obj = PlotResidualPTps(edi_file_list_1, edi_file_list_2)
        plt_obj.plot()
        plt.pause(1)
