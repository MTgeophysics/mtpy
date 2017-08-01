import glob
import inspect
import unittest
from unittest import TestCase

# configure matplotlib for testing
import matplotlib.pyplot as plt

plt.ion()
import os

from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection

edi_paths = [
    "",
    "tests\\data\\edifiles",
    "examples\\data\\edi2",
    "examples\\data\\edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests\\data\\edifiles2"
]


class TestPlotResPhasePseudoSection(TestCase):
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
        edi_path = edi_paths[1]
        self._plot(edi_path, "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_plot_02(self):
        edi_path = edi_paths[2]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03(self):
        edi_path = edi_paths[3]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    def test_plot_04(self):
        edi_path = edi_paths[4]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    def test_plot_05(self):
        edi_path = edi_paths[5]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    def _plot(self, edi_path, save_figure_path):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        save_figure_path = os.path.join(self._temp_dir, save_figure_path)
        pt_obj = PlotResPhasePseudoSection(fn_list=edi_file_list, plot_yn='n', plot_style='imshow')
        pt_obj.plot()
        # pt_obj.save_plot(save_figure_path)
        plt.pause(1)
        pt_obj = PlotResPhasePseudoSection(fn_list=edi_file_list, plot_yn='n', plot_style='pcolormesh')
        pt_obj.plot()
        # pt_obj.save_plot(save_figure_path)
        plt.pause(1)
