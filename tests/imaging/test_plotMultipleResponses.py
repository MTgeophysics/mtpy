import glob
import os
import unittest
from unittest import TestCase

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update(matplotlib.rcParamsDefault)  # reset matplotlib params
plt.ion()

from mtpy.imaging.plotnresponses import PlotMultipleResponses
from mtpy.utils.decorator import ImageCompare


edi_paths = [
    "",
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests/data/edifiles2"
]


class TestPlotMultipleResponses(TestCase):
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

    @ImageCompare(fig_size=(16, 6),  savefig_kwargs={'dpi': 200})
    def test_plot_01_all(self):
        edi_path = edi_paths[1]
        self._plot_all(edi_path)

    @ImageCompare(fig_size=(8, 6), savefig_kwargs={'dpi': 200})
    def test_plot_01_compare(self):
        edi_path = edi_paths[1]
        self._plot_compare(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    @ImageCompare(fig_size=(16, 6), savefig_kwargs={'dpi': 200})
    def test_plot_02_all(self):
        edi_path = edi_paths[2]
        self._plot_all(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    @ImageCompare(fig_size=(8, 6), savefig_kwargs={'dpi': 200})
    def test_plot_02_compare(self):
        edi_path = edi_paths[2]
        self._plot_compare(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    @ImageCompare(fig_size=(16, 6), savefig_kwargs={'dpi': 200})
    def test_plot_03_all(self):
        edi_path = edi_paths[3]
        self._plot_all(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    @ImageCompare(fig_size=(8, 6), savefig_kwargs={'dpi': 200})
    def test_plot_03_compare(self):
        edi_path = edi_paths[3]
        self._plot_compare(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    @ImageCompare(fig_size=(16, 6), savefig_kwargs={'dpi': 200})
    def test_plot_04_all(self):
        edi_path = edi_paths[4]
        self._plot_all(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    @ImageCompare(fig_size=(8, 6), savefig_kwargs={'dpi': 200})
    def test_plot_04_compare(self):
        edi_path = edi_paths[4]
        self._plot_compare(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    @ImageCompare(fig_size=(32, 6), savefig_kwargs={'dpi': 200})
    def test_plot_05_all(self):
        edi_path = edi_paths[5]
        self._plot_all(edi_path)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    @ImageCompare(fig_size=(8, 6), savefig_kwargs={'dpi': 200})
    def test_plot_05_compare(self):
        edi_path = edi_paths[5]
        self._plot_compare(edi_path)

    @staticmethod
    def _plot_all(edi_path):
        edi_file_list = glob.glob(os.path.join(edi_path, '*.edi'))
        # 1
        # for file in edi_file_list:
        #     if self._progress_bar:
        #         self._progress_bar.onStart()
        #     pt_obj = PlotMultipleResponses(
        #         fn_list=[file],
        #         plot_num=1,
        #         plot_tipper='yr',
        #         plot_style='1'
        #     )
        #     pt_obj.plot()
        #     if self._progress_bar:
        #         self._progress_bar.onFinished()
        #     plt.pause(0.5)
        # all
        pt_obj = PlotMultipleResponses(
            fn_list=edi_file_list,
            plot_num=1,
            plot_tipper='yr',
            plot_style='all',
            plot_yn='n'
        )
        pt_obj.plot()
        plt.pause(0.5)

    @staticmethod
    def _plot_compare(edi_path):
        edi_file_list = glob.glob(os.path.join(edi_path, '*.edi'))
        # compare
        pt_obj = PlotMultipleResponses(
            fn_list=edi_file_list,
            plot_num=1,
            plot_tipper='yr',
            plot_style='compare',
            plot_yn='n'
        )
        pt_obj.plot()
        plt.pause(0.5)
