
import glob
import inspect
import os
import unittest
from unittest import TestCase

# configure matplotlib for testing
import matplotlib.pyplot as plt

from mtpy.imaging.plotstrike2d import PlotStrike2D

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


class TestPlotStrike2D(TestCase):
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
        pt_obj = PlotStrike2D(fn_list=edi_file_list, plot_yn='n')
        pt_obj.plot()
        plt.pause(1)

        # change rotation
        pt_obj = PlotStrike2D(fn_list=edi_file_list, plot_yn='n', rot_z=90)
        pt_obj.plot()
        plt.pause(1)

        # plot type
        pt_obj = PlotStrike2D(fn_list=edi_file_list, plot_yn='n', plot_type=1)
        pt_obj.plot()
        plt.pause(1)

        # plot_tipper
        pt_obj = PlotStrike2D(fn_list=edi_file_list, plot_yn='n', plot_tipper='y')
        pt_obj.plot()
        plt.pause(1)

        # fold
        pt_obj = PlotStrike2D(fn_list=edi_file_list, plot_yn='n', fold=False)
        pt_obj.plot()
        plt.pause(1)
